"""Deep Potential model / ABACUS-DP capability detection.

Before generating a DP-MD INPUT two *independent* questions must be answered:

1. What can this ABACUS binary read?  — a compile-time property.
2. What format is this model file?    — a property of the file alone.

For (1) the conda env on PATH is **not** the answer.  ABACUS links the
DeepMD-kit C++ API at build time (e.g. against ``deepmd-2.2.11-cxx``), so
``conda activate deepmd_v3`` does not widen what an existing binary can
read.  Asking the binary itself is the only ground truth that does not
depend on the user having configured anything correctly.

Probing is layered, most trustworthy first:

* ``explicit``     — ``paths.abacus_dp_version`` in config (a human said so)
* ``cmake_cache``  — ``build/CMakeCache.txt`` of the source tree, if kept
* ``binary_probe`` — ``strings`` + ``ldd`` on the binary itself (the default)
* ``env_hint``     — ``paths.slurm_env_file``; a *declaration*, not a fact.
                      Recorded for cross-checking only, never as the version.

Results are cached per host in ``~/.abacuscopilot/machines.yaml`` because the
answer only changes when someone rebuilds ABACUS.
"""

from __future__ import annotations

import re
import shutil
import socket
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# Markers
# ---------------------------------------------------------------------------

#: Node carrying the type map of a DeepMD-kit v2 TensorFlow frozen graph.
_TMAP_MARKER = b"model_attr/tmap"
#: Node carrying the frozen-graph format version (e.g. b"1.1").
_VERSION_NODE = b"model_attr/model_version"
#: Any of these means the bytes are a DeepMD frozen TensorFlow graph.  Needed
#: separately from ``_TMAP_MARKER`` so "a TF graph whose tmap was lost" can be
#: told apart from "not a TF graph at all" — the former is the signature of a
#: model that a converter has mangled, and is worth naming precisely.
_TF_GRAPH_MARKERS = (b"descrpt_attr/", b"fitting_attr/", b"model_attr/model_type")

#: Symbols / TF op names that only appear when ABACUS really links DeePMD.
#: A bare "deepmd" is deliberately NOT a positive marker: ABACUS also ships a
#: "please recompile with DeePMD" style message, so the word alone proves
#: nothing about the build.
_DP_POSITIVE_MARKERS = (
    "deepmd::DeepPot",
    "DeepPotTF",
    "ProdEnvMatA",
    "ProdForceSeA",
    "libdeepmd",
)
#: Co-occurring with "deepmd", these indicate a binary built WITHOUT DP.
_DP_ABSENT_MARKERS = (
    "USE_DEEPMD",
    "not compiled",
    "recompile",
    "not enabled",
)

_VERSION_RE = re.compile(r"(\d+\.\d+\.\d+)")

#: A version only counts when it sits next to the deepmd name itself.
#:
#: A bare X.Y.Z search over a ``strings`` dump of a debug-info binary finds
#: build paths, not versions.  ABACUS's own binary embeds, on one line,
#: ``/home/young/miniconda3/envs/deepmd_v2/lib:/…/abacus-develop-LTSv3.10.0/…``
#: — the word "deepmd" is in there (the env name), and the first X.Y.Z on the
#: line is 3.10.0, ABACUS's version.  The probe would then report DeepMD-kit
#: as 3.10.0 and, because the compatibility verdict branches on the major
#: version, warn about a perfectly healthy setup.
#:
#: ``[^0-9]{0,20}`` keeps the usual spellings working — ``deepmd-2.2.11``,
#: ``libdeepmd.so.2.2.11``, ``DeePMD-kit version: 2.2.11`` — while refusing to
#: skip over a digit, so a path like ``envs/deepmd_v2/lib`` can never reach a
#: version further along the line.
_DEEPMD_VERSION_RE = re.compile(r"deepmd[^0-9]{0,20}(\d+\.\d+\.\d+)", re.IGNORECASE)

#: Suffixes ABACUS's ``dp`` esolver can load.  ``.pb`` is the TensorFlow graph
#: (DeepMD v1/v2), the rest are PyTorch-format models.
_TORCH_SUFFIXES = (".pt", ".pth", ".dpmodel")
MODEL_SUFFIXES = (".pb", *_TORCH_SUFFIXES)

#: Statuses returned by :func:`assess_compatibility`.
OK = "ok"
WARN = "warn"
NEEDS_CONVERSION = "needs_conversion"
INCOMPATIBLE = "incompatible"
UNKNOWN = "unknown"
MISSING = "missing"


# ---------------------------------------------------------------------------
# Model side
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ModelFingerprint:
    """What a potential file actually is, determined from its bytes."""

    path: str
    exists: bool
    kind: str  # tf_graphdef | torch_model | dpmodel_zip | empty | missing | unknown
    has_tmap: bool = False
    model_version: str | None = None
    size: int = 0

    @property
    def is_dp2_tf(self) -> bool:
        """A TensorFlow graph carrying the descriptors ABACUS reads."""
        return self.kind == "tf_graphdef" and self.has_tmap

    @property
    def has_version(self) -> bool:
        """Whether ``model_attr/model_version`` is present.

        Not cosmetic. DeePMD-kit's C API reads a missing version node as
        ``0.0`` and refuses anything below ``1.1``, so a graph that looks
        perfectly healthy — type map and all — aborts the run with
        "incompatable model: version 0.0 in graph, but version 1.1 supported".
        Models published in the older, uncompressed DP format look exactly
        like this; ``dp compress`` is what stamps the node.
        """
        return self.model_version is not None

    @property
    def is_torch(self) -> bool:
        """A PyTorch-format model — readable only by a PyTorch-backend build.

        Same zip container as ``dpmodel_zip``; the file suffix is what tells a
        DP3 ``.pt`` model apart from a ``.pb`` that is really an archive.
        """
        return self.kind == "torch_model"

    @property
    def usable(self) -> bool:
        """Whether a DeePMD ABACUS build could plausibly load this.

        Used to decide whether offering a backup model is even worthwhile —
        distinct from :meth:`Compatibility.is_problem`, which also covers
        uncertainty about the binary.
        """
        return self.is_dp2_tf and self.has_version

    def describe(self) -> str:
        if not self.exists:
            return f"{self.path}: 文件不存在"
        if self.kind == "empty":
            return f"{self.path}: 空文件 (0 字节)"
        if self.kind == "dpmodel_zip":
            return (f"{self.path}: DeepMD v3 后端格式 (zip), "
                    f"{self.size} 字节")
        if self.kind == "torch_model":
            return (f"{self.path}: PyTorch 格式模型 (zip), "
                    f"{self.size} 字节")
        if self.kind == "tf_graphdef":
            tmap = "有 type map" if self.has_tmap else "缺少 model_attr/tmap"
            ver = (f", model_version={self.model_version}" if self.has_version
                   else ", 无 model_version 节点 —— ABACUS 会拒收")
            return f"{self.path}: TensorFlow GraphDef ({tmap}{ver}), {self.size} 字节"
        return f"{self.path}: 无法识别的格式, {self.size} 字节"


def discover_models(directory: str | Path = ".") -> list[ModelFingerprint]:
    """Every DP model file sitting directly in *directory*, best candidate first.

    The template default (``graph-compress.pb``) is only a default — most
    people keep their model under whatever name the training run produced, and
    renaming it to suit the tool is a step that shouldn't exist.  So the input
    generator looks for itself and asks.

    Ordered by "can ABACUS probably load this" first, then by name, so the head
    of the list is a usable default.  Broken models are still returned: when
    the directory's only ``.pb`` is the damaged one, dropping it would leave the
    user looking at an empty list instead of the file they need to repair.
    """
    root = Path(directory).expanduser()
    if not root.is_dir():
        return []

    paths = sorted(
        (p for p in root.iterdir()
         if p.is_file() and p.suffix.lower() in MODEL_SUFFIXES),
        key=lambda p: p.name,
    )
    prints = [fingerprint_model(p) for p in paths]
    prints.sort(key=lambda fp: (not fp.usable, fp.path))
    return prints


def _scan_markers(path: Path, markers: tuple[bytes, ...], *,
                  chunk_size: int = 1 << 20, overlap: int = 64) -> dict[bytes, int]:
    """One streaming pass; maps each found marker to its first file offset.

    Models run to hundreds of MB, so this never holds the file in memory and
    stops as soon as every marker has been seen.
    """
    wanted = set(markers)
    found: dict[bytes, int] = {}
    tail = b""
    base = 0
    with open(path, "rb") as fh:
        while wanted:
            block = fh.read(chunk_size)
            if not block:
                break
            buf = tail + block
            for marker in list(wanted):
                idx = buf.find(marker)
                if idx >= 0:
                    found[marker] = base - len(tail) + idx
                    wanted.discard(marker)
            tail = buf[-overlap:]
            base += len(block)
    return found


def _extract_model_version(path: Path, offset: int | None) -> str | None:
    """Best-effort read of the frozen-graph version.

    The value sits in the (uncompressed) GraphDef bytes shortly after the node
    name.  Protobuf length prefixes are not ASCII digits, so a plain digit-run
    regex over a tight window is reliable enough — this is a diagnostic aid,
    not a correctness-critical value.
    """
    if offset is None:
        return None
    try:
        with open(path, "rb") as fh:
            fh.seek(offset)
            window = fh.read(160)
    except OSError:
        return None
    match = re.search(rb"([0-9]+\.[0-9]+(?:\.[0-9]+)?)", window)
    return match.group(1).decode() if match else None


def fingerprint_model(path: str | Path) -> ModelFingerprint:
    """Identify a potential file without importing deepmd or tensorflow.

    Purely structural: a ``.pb`` is either a TensorFlow GraphDef (contains
    ``model_attr/tmap``) or a DeepMD v3 backend archive (a zip, ``PK``).
    """
    p = Path(path).expanduser()
    if not p.is_file():
        return ModelFingerprint(path=str(p), exists=False, kind="missing")

    size = p.stat().st_size
    if size == 0:
        return ModelFingerprint(path=str(p), exists=True, kind="empty", size=0)

    try:
        with open(p, "rb") as fh:
            magic = fh.read(2)
    except OSError:
        return ModelFingerprint(path=str(p), exists=True, kind="unknown", size=size)

    if magic == b"PK":
        # A zip, but which kind depends on the name: DP3's PyTorch checkpoints
        # are zips too, and calling one a "v3 backend archive" would send the
        # user off to convert a file that was already in the right format.
        kind = "torch_model" if p.suffix.lower() in _TORCH_SUFFIXES else "dpmodel_zip"
        return ModelFingerprint(path=str(p), exists=True, kind=kind, size=size)

    found = _scan_markers(p, (*_TF_GRAPH_MARKERS, _TMAP_MARKER, _VERSION_NODE))
    is_tf_graph = any(m in found for m in _TF_GRAPH_MARKERS)
    return ModelFingerprint(
        path=str(p),
        exists=True,
        kind="tf_graphdef" if is_tf_graph else "unknown",
        has_tmap=_TMAP_MARKER in found,
        model_version=_extract_model_version(p, found.get(_VERSION_NODE)),
        size=size,
    )


# ---------------------------------------------------------------------------
# Binary side
# ---------------------------------------------------------------------------


@dataclass
class DPCapability:
    """What an ABACUS binary was compiled to do, and how we know."""

    binary: str
    resolved: str | None = None
    has_dp: bool | None = None          # None = could not determine
    deepmd_version: str | None = None
    linkage: str = "unknown"            # static | dynamic | unknown
    linked_libs: list[str] = field(default_factory=list)
    runtime_env: str = ""               # hint from config, not a fact
    source: str = "unknown"             # provenance of deepmd_version/has_dp
    probed_at: str = ""

    def describe(self) -> str:
        if self.resolved is None:
            return f"{self.binary}: 未找到该二进制"
        if self.has_dp is True:
            ver = f" DeepMD {self.deepmd_version}" if self.deepmd_version else " (版本未知)"
            return f"{self.resolved}: 支持 DP{ver} [{self.linkage}, 来源 {self.source}]"
        if self.has_dp is False:
            return f"{self.resolved}: 未编译 DP 支持"
        return f"{self.resolved}: DP 支持情况未知 (来源 {self.source})"


def _run(cmd: list[str], timeout: int = 30) -> str:
    """Run *cmd* and return stdout+stderr, or "" when unavailable."""
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return (proc.stdout or "") + (proc.stderr or "")
    except (OSError, subprocess.SubprocessError):
        return ""


def _classify_dp_support(blob: str) -> bool | None:
    """True / False / None from a *strings* dump of an ABACUS binary."""
    low = blob.lower()
    if any(marker.lower() in low for marker in _DP_POSITIVE_MARKERS):
        return True
    if "deepmd" in low and any(m.lower() in low for m in _DP_ABSENT_MARKERS):
        return False
    return None


def _parse_version(lines: list[str]) -> str | None:
    """DeepMD version from *lines*, or None.

    Only counts a version welded to the deepmd name — see
    :data:`_DEEPMD_VERSION_RE` for why a bare ``X.Y.Z`` search is not enough.
    """
    for line in lines:
        match = _DEEPMD_VERSION_RE.search(line)
        if match:
            return match.group(1)
    return None


def _version_from_linked_lib(libs: list[str]) -> str | None:
    """DeepMD version read from the env that owns a linked ``libdeepmd``.

    ``libdeepmd_c.so`` carries no version string of its own, so scanning it is
    a dead end.  The conda env around it is not: a library at
    ``<prefix>/lib/libdeepmd_c.so`` implies ``<prefix>`` also holds
    ``site-packages/deepmd_kit-<version>.dist-info``, whose METADATA states the
    version outright.  Pure file reads — no subprocess, no import.
    """
    for lib in libs:
        path = Path(lib)
        if "deepmd" not in path.name.lower():
            continue
        prefix = path.parent.parent if path.parent.name in ("lib", "lib64") else path.parent
        for meta in prefix.glob("lib*/python*/site-packages/deepmd_kit-*.dist-info/METADATA"):
            for line in meta.read_text(errors="replace").splitlines():
                if line.startswith("Version:"):
                    return line.split(":", 1)[1].strip()
    return None


def _linkage_info(lines: list[str]) -> tuple[str, list[str]]:
    """(linkage, matching lib paths) from ``ldd`` / ``otool -L`` output."""
    libs = [ln.strip() for ln in lines
            if re.search(r"deepmd|tensorflow", ln, re.IGNORECASE)]
    if not libs:
        return "unknown", []
    # "libdeepmd.so => /path/to/libdeepmd.so (0x...)" -> keep the real path
    paths = []
    for line in libs:
        for token in line.split():
            if token.startswith("/") and "=>" not in token:
                paths.append(token)
                break
        else:
            paths.append(line)
    return "dynamic", paths


def _resolve_binary(binary: str | None, config: dict) -> str:
    if binary:
        return binary
    paths = config.get("paths", {}) if isinstance(config, dict) else {}
    return (paths.get("abacus_dp_binary")
            or paths.get("abacus_binary")
            or "abacus")


def _probe_cmake_cache(config: dict) -> tuple[bool | None, str | None]:
    """(has_dp, version) from the source tree's CMakeCache, if it exists."""
    src = (config.get("paths", {}) or {}).get("abacus_source") or ""
    if not src:
        return None, None
    root = Path(src).expanduser()
    if not root.is_dir():
        return None, None
    for cache in list(root.glob("build*/CMakeCache.txt")) + list(root.glob("CMakeCache.txt")):
        text = _run(["grep", "-E", "USE_DEEPMD|DEEPMD", str(cache)], timeout=15)
        if not text:
            continue
        on = re.search(r"USE_DEEPMD[^=]*=\s*(ON|TRUE|1)", text, re.IGNORECASE)
        if not on:
            continue
        version = _parse_version(text.splitlines())
        return True, version
    return None, None


def probe_abacus_dp(binary: str | None = None, *, config: dict | None = None,
                    use_cache: bool = True, refresh: bool = False) -> DPCapability:
    """Determine whether *binary* can do DP, and against which DeepMD version.

    Reads the host cache unless *refresh* is set.  Never raises: an
    undeterminable capability comes back with ``has_dp=None``.
    """
    if config is None:
        from abacuscopilot.config import load_config
        config = load_config()

    name = _resolve_binary(binary, config)
    resolved = shutil.which(name)
    if resolved is None:
        candidate = Path(name).expanduser()
        resolved = str(candidate) if candidate.is_file() else None

    if resolved is None:
        return DPCapability(binary=name, source="not_found",
                            probed_at=date.today().isoformat())

    if use_cache and not refresh:
        cached = _load_cached_capability(resolved)
        if cached is not None:
            return cached

    cap = DPCapability(binary=name, resolved=resolved,
                       probed_at=date.today().isoformat())

    # --- L1: the build cache is the compile-time record itself -------------
    has_dp, version = _probe_cmake_cache(config)
    if has_dp is not None:
        cap.has_dp, cap.deepmd_version, cap.source = has_dp, version, "cmake_cache"

    # --- L2: ask the binary ----------------------------------------------
    if cap.has_dp is None or cap.deepmd_version is None:
        blob = _run(["strings", resolved], timeout=60)
        if cap.has_dp is None:
            cap.has_dp = _classify_dp_support(blob)
            if cap.has_dp is not None:
                cap.source = "binary_probe"
        if cap.deepmd_version is None:
            cap.deepmd_version = _parse_version(blob.splitlines())

        linkage_out = _run(["ldd", resolved], timeout=30) or _run(
            ["otool", "-L", resolved], timeout=30)
        if linkage_out:
            cap.linkage, cap.linked_libs = _linkage_info(linkage_out.splitlines())
            if cap.linkage == "dynamic":
                cap.deepmd_version = cap.deepmd_version or _parse_version(
                    cap.linked_libs) or _version_from_linked_lib(cap.linked_libs)
        if cap.linkage == "unknown" and cap.has_dp is True:
            cap.linkage = "static"

    # --- L4 / L0: declarations (hint, then human override) ----------------
    cap.runtime_env = (config.get("paths", {}) or {}).get("slurm_env_file", "") or ""
    override = ((config.get("paths", {}) or {}).get("abacus_dp_version") or "").strip()
    if override:
        cap.deepmd_version = override
        cap.source = "explicit" if cap.source in ("unknown", "not_found") \
            else f"{cap.source}+explicit"

    save_machine_capability(cap)
    return cap


def check_runtime_env_hint(cap: DPCapability) -> str | None:
    """Warn when a dynamically linked binary has no env script recorded.

    A dynamic ABACUS needs the matching ``libdeepmd`` / ``libtensorflow`` on
    ``LD_LIBRARY_PATH`` at run time; ``paths.slurm_env_file`` is where the
    user declares that.
    """
    if cap.linkage == "dynamic" and not cap.runtime_env:
        return ("该 abacus 动态链接 DeepMD/TensorFlow,但 config 里没配 "
                "paths.slurm_env_file —— 作业里可能找不到库。")
    return None


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Compatibility:
    status: str
    message: str
    guidance: str = ""

    @property
    def is_problem(self) -> bool:
        """True unless this is a clean pass — WARN included.

        A warn still means "not confirmed good", so callers must not render it
        as a green check.
        """
        return self.status != OK


def assess_compatibility(model: ModelFingerprint,
                         cap: DPCapability) -> Compatibility:
    """Pair a model fingerprint with a binary capability."""
    if not model.exists:
        return Compatibility(MISSING, "模型文件不存在。",
                             "把 DeePMD 模型放到当前目录,或改 INPUT 的 pot_file。")

    if model.kind == "empty":
        return Compatibility(
            INCOMPATIBLE, "模型文件是 0 字节 —— 很可能是一次失败的转换留下的残片。",
            "换回可用的模型文件;若旁边有 *_original 备份,优先试它。")

    if cap.has_dp is False:
        return Compatibility(
            INCOMPATIBLE, "该 abacus 没有编译 DeepMD 支持,无法做 DP 计算。",
            "改用编译了 DP 的二进制(如 abacus-LTSv3.10.0-dp),或在 config 的 "
            "paths.abacus_dp_binary 指向它。")

    if model.kind == "dpmodel_zip":
        return Compatibility(
            NEEDS_CONVERSION, "模型是 DeepMD v3 后端格式 (zip),该 abacus 读不了。",
            "用 dp convert-backend 转成 TF .pb,或换用支持该格式的二进制。")

    if model.kind == "torch_model":
        # PyTorch models need a PyTorch-backend ABACUS.  The symbol probe only
        # answers "is DP linked at all", not which backend, so this stays a
        # warning: converting a working .pt would be worse than trying it.
        return Compatibility(
            WARN, "模型是 PyTorch 格式 (.pt),只有 PyTorch 后端的 abacus 能读。",
            "若该 abacus 编译时用的是 TF 后端,需要 dp convert-backend 转成 .pb;"
            "先跑一步试算确认。")

    if model.kind == "tf_graphdef" and not model.has_tmap:
        return Compatibility(
            INCOMPATIBLE,
            "TensorFlow 图里没有 model_attr/tmap 节点 —— ABACUS 会读到空 type map,"
            "然后对每个元素报 \"not found in the type map\"。",
            "这个文件已经被某个转换工具改坏了。换回可用的模型;若旁边有 "
            "*_original 备份,优先试它。")

    if model.kind == "tf_graphdef" and not model.has_version:
        return Compatibility(
            INCOMPATIBLE,
            "模型缺 model_attr/model_version 节点 —— DeePMD-kit 按 version 0.0 处理,"
            "会直接 abort:\"incompatable model: version 0.0 in graph, but version 1.1 "
            "supported\"。type map 齐全也照样跑不起来。",
            "这是旧格式(未压缩)的模型。用 dp compress 重新压缩一遍,它会给图打上 "
            "model_version 节点;例如:dp compress -i 原模型.pb -o graph-compress.pb")

    if not model.is_dp2_tf:
        return Compatibility(UNKNOWN, "无法识别的模型格式。", model.describe())

    # A proper DeepMD v2 TF graph.  Whether the binary likes it depends on the
    # DeepMD version it was built against.
    major = (cap.deepmd_version or "").split(".")[0]
    if cap.has_dp is None:
        return Compatibility(
            WARN, "模型是 DeepMD v2 TF 图,但没能确定 abacus 的 DP 能力。",
            "跑一次试算确认;或在 config 里显式设置 paths.abacus_dp_version。")
    if major == "2" or not major:
        return Compatibility(
            OK, f"模型格式匹配 (DeepMD v2 TF 图{', 版本 ' + cap.deepmd_version if cap.deepmd_version else ''})。")
    return Compatibility(
        WARN, f"模型是 DeepMD v2 TF 图,而 abacus 编译时用的是 DeepMD {cap.deepmd_version}。",
        "多数情况下仍可运行(DP3 保留了 TF 后端),建议先跑一步试算确认。")


# ---------------------------------------------------------------------------
# Per-host cache
# ---------------------------------------------------------------------------


def _machines_path() -> Path:
    return Path.home() / ".abacuscopilot" / "machines.yaml"


def _load_cached_capability(resolved_binary: str) -> DPCapability | None:
    path = _machines_path()
    if not path.is_file():
        return None
    try:
        data = yaml.safe_load(path.read_text()) or {}
    except (OSError, yaml.YAMLError):
        return None
    entry = (data.get("hosts", {}) or {}).get(socket.gethostname(), {}) or {}
    stored = entry.get("abacus_dp")
    if not isinstance(stored, dict) or stored.get("resolved") != resolved_binary:
        return None
    known = {f for f in DPCapability.__dataclass_fields__}
    return DPCapability(**{k: v for k, v in stored.items() if k in known})


def save_machine_capability(cap: DPCapability) -> None:
    """Record *cap* under this host.  Failures are non-fatal by design."""
    path = _machines_path()
    try:
        data = yaml.safe_load(path.read_text()) if path.is_file() else {}
        if not isinstance(data, dict):
            data = {}
        data.setdefault("hosts", {})
        data["hosts"].setdefault(socket.gethostname(), {})
        data["hosts"][socket.gethostname()]["abacus_dp"] = asdict(cap)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True))
    except (OSError, yaml.YAMLError):
        pass


def get_abacus_dp_capability(*, binary: str | None = None, config: dict | None = None,
                            refresh: bool = False) -> DPCapability:
    """Cached capability for this host (probes once, then reuses)."""
    return probe_abacus_dp(binary, config=config, use_cache=True, refresh=refresh)
