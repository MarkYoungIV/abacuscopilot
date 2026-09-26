"""Configuration management for abacuscopilot.

Reads and writes YAML config from ~/.abacuscopilot/config.yaml.
Settings include default paths, pseudopotential directories,
plotting preferences, and user-defined presets.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def release_root() -> Path:
    """The release tree holding ``pyproject.toml``, ``PP-Orb/`` and ``scripts/``.

    Found by walking up to the ``pyproject.toml`` marker rather than by
    counting ``.parent`` hops from ``__file__``.  Hop counting is what silently
    pointed ``PP-Orb/`` at ``<release>/src`` the moment the package moved under
    ``src/`` — the count depends on the layout, and every copy of it has to be
    found and fixed separately.  The marker does not.

    A non-editable install (site-packages, no release tree above it) has no such
    marker; the package's parent directory is returned, where detection finds no
    ``PP-Orb/`` and correctly yields nothing.
    """
    pkg_dir = Path(__file__).resolve().parent
    # Three levels up covers both layouts (<release>/abacuscopilot and
    # <release>/src/abacuscopilot); the bound keeps a stray marker in $HOME
    # from being picked up.
    for candidate in (pkg_dir, *list(pkg_dir.parents)[:3]):
        if (candidate / "pyproject.toml").is_file():
            return candidate
    return pkg_dir.parent


def _path_within(path: Path, root: Path) -> bool:
    """True when *path* equals *root* or lives somewhere under it."""
    try:
        path.expanduser().resolve().relative_to(root.expanduser().resolve())
        return True
    except ValueError:
        return False


def _registered_family_dirs(config: dict) -> list[Path]:
    """Resolved dirs the user registered as external-series roots.

    Walks ``libraries.families.<id>.pseudo_dir`` plus every
    ``...orbital_dirs.*`` entry.  Used to keep a registered family (e.g. an
    ABACUS-APNS copy dropped under ``PP-Orb/``) out of the auto-detected
    "bundled" SG15 lists, so the two series never silently mix.
    """
    out: list[Path] = []
    libs = config.get("libraries", {})
    fams = libs.get("families", {})
    if not isinstance(fams, dict):
        return out
    for fam in fams.values():
        if not isinstance(fam, dict):
            continue
        for v in (fam.get("pseudo_dir"),):
            if isinstance(v, str) and v:
                out.append(Path(v).expanduser())
        orb = fam.get("orbital_dirs")
        if isinstance(orb, dict):
            for v in orb.values():
                if isinstance(v, str) and v:
                    out.append(Path(v).expanduser())
    return out


# ---------------------------------------------------------------------------
# The PP-Orb/ layout, as declared data
# ---------------------------------------------------------------------------
#
# These two tables are the single source of truth for what the auto-detection
# is allowed to see.  Nothing infers a series from the shape of the tree: the
# previous rule ("every top-level folder under PP-Orb/ is a library root") read
# one folder per *library*, which stopped being true once a series became one
# folder holding both libraries — and, worse, it admitted the external series
# folders that physically live under PP-Orb/ into the bundled default lists, so
# an SG15 job silently got Dojo pseudopotentials.
#
# One top-level folder per series.  ``id -> (top folder, {extension: the
# sub-path under it holding that extension's files})``.  A series whose .upf
# and .orb sit in the same folder repeats the path.
_BUNDLED_SERIES_LAYOUT: dict[str, tuple[str, dict[str, str]]] = {
    "sg15": (
        "SG15-Version1p0",
        {
            ".upf": "SG15-Version1p0_Pseudopotential",
            ".orb": "SG15-Version1p0__StandardOrbitals-Version2p0",
        },
    ),
    "lanthanides": (
        "lanthanides-f--core.icmod1",
        {
            # Large-core (f-electron-pseudized) 4f PPs and their NAOs sit side
            # by side under PD04.3+f--core.icmod1/{El}/, so both extensions
            # share one root.  Zero element overlap with SG15 (14 elements,
            # Ce..Lu, none of which SG15's orbitals cover), so listing this
            # series alongside SG15 can never shadow an SG15 element — it only
            # fills the 4f gap, and La stays pseudopotential-only.
            ".upf": "PD04.3+f--core.icmod1",
            ".orb": "PD04.3+f--core.icmod1",
        },
    ),
}

# External series (NOT bundled): the top folder a user's own copy of the series
# lives under when dropped into PP-Orb/, plus the sub-paths to its pseudo and
# (per variant) orbital roots.  Used to auto-register a series that is already
# physically present instead of prompting for its paths by hand — and to keep
# such a folder OUT of the bundled lists even before it is registered, which is
# the bug this table was extended for.  A shared orbital root (Dojo keeps
# SZ/DZP/TZDP all under one Orbitals_v2.0) repeats the same sub-path.
_EXTERNAL_SERIES_LAYOUT: dict[str, tuple[str, str, dict[str, str]]] = {
    "apns": (
        "ABACUS-APNS-PPORBs-v1",
        "apns-pseudopotentials-v1",
        {
            "efficiency": "apns-orbitals-efficiency-v1",
            "precision": "apns-orbitals-precision-v1",
        },
    ),
    "dojoncfr": (
        "Dojo-NC-FR",
        "Pseudopotential",
        {"sz": "Orbitals_v2.0", "dzp": "Orbitals_v2.0", "tzdp": "Orbitals_v2.0"},
    ),
}


def _bundled_series_root(pp_orb: Path, series: str, ext: str) -> Path | None:
    """The sub-root one bundled *series* declares for *ext*, if it exists.

    ``None`` when the series declares no path for that extension (e.g. the
    lanthanide supplement has no separate orbital root) or when the folder is
    absent.  There is deliberately no fallback to the series' own folder: it
    holds both extensions, so falling back would collapse ``pseudo_library``
    and ``orbital_library`` onto one root — the bug this replaced.
    """
    top, subs = _BUNDLED_SERIES_LAYOUT[series]
    sub = subs.get(ext)
    if sub is None:
        return None
    d = pp_orb / top / sub
    return d.resolve() if d.is_dir() else None


def _bundled_series_roots(pp_orb: Path, ext: str) -> list[Path]:
    """Existing sub-roots of every bundled series under *pp_orb* holding *ext*.

    One entry per series that both declares a path for *ext* and has it on disk.
    """
    return [
        root
        for series in _BUNDLED_SERIES_LAYOUT
        if (root := _bundled_series_root(pp_orb, series, ext)) is not None
    ]


def _bundled_series_tops() -> set[str]:
    """Top-level ``PP-Orb/`` folders owned by a bundled series."""
    return {top for top, _ in _BUNDLED_SERIES_LAYOUT.values()}


def _external_series_tops() -> set[str]:
    """Top-level ``PP-Orb/`` folders owned by a known external series.

    Derived from :data:`_EXTERNAL_SERIES_LAYOUT` so that a known series folder
    is kept out of the bundled lists immediately — registration under
    ``libraries.families`` is *not* required, and on a machine that never
    registered it the folder must still not be searched as SG15.
    """
    return {top for top, _, _ in _EXTERNAL_SERIES_LAYOUT.values()}


def _detect_library_dirs(ext: str, exclude_containing: object = ()) -> list[str]:
    """Auto-detect the bundled library directories inside the package root.

    The bundled series are declared in :data:`_BUNDLED_SERIES_LAYOUT`, and each
    one contributes the explicit sub-root that holds *ext* — never its series
    folder, which holds both extensions at once::

      PP-Orb/SG15-Version1p0/SG15-Version1p0_Pseudopotential/              (.upf)
      PP-Orb/SG15-Version1p0/SG15-Version1p0__StandardOrbitals-Version2p0/ (.orb)
      PP-Orb/lanthanides-f--core.icmod1/PD04.3+f--core.icmod1/             (both)

    The legacy ``Pseudopotential/`` / ``Orbitals/`` pair is still honoured when
    there is no ``PP-Orb/`` directory at all.

    Two kinds of top-level folder are deliberately *not* searched: a bundled
    series folder (represented by the sub-roots above), and a known external
    series folder (``ABACUS-APNS-PPORBs-v1``, ``Dojo-NC-FR``) — those series
    must never leak into the default SG15 lists, whether or not they are
    registered under ``libraries.families``.  Anything else — a folder a user
    dropped in themselves — is still detected when it holds *ext anywhere below
    it, which is the drop-in behaviour users are told to rely on.

    *exclude_containing*: iterable of dirs (typically from
    :func:`_registered_family_dirs`).  A candidate root that *contains* (or
    equals) one of them is skipped, so a registered family is never absorbed
    into the bundled lists even when its top folder is not a known one.  The
    check is applied uniformly, declared roots included: narrowing it would
    re-admit a foreign series into the bundled list for a user who nests one
    there on purpose.

    Returns a list of absolute paths (possibly empty).
    """
    pkg_root = release_root()
    pp_orb = pkg_root / "PP-Orb"
    if pp_orb.is_dir():
        skip = _bundled_series_tops() | _external_series_tops()
        dropins = [d for d in sorted(pp_orb.iterdir())
                   if d.is_dir() and not d.name.startswith(".") and d.name not in skip]
        roots = _bundled_series_roots(pp_orb, ext) + dropins
    else:
        roots = [d for name in ("Pseudopotential", "Orbitals")
                 if (d := pkg_root / name).is_dir()]

    exclude = list(exclude_containing or ())

    found: list[str] = []
    for root in roots:
        if any(_path_within(p, root) for p in exclude):
            continue  # a registered family lives here — don't bundle it
        if any(f.is_file() and f.name.lower().endswith(ext)
               for f in root.rglob("*")):
            found.append(str(root.resolve()))
    return found


DEFAULT_CONFIG = {
    "defaults": {
        "pseudo_dir": "./",
        "orbital_dir": "./",
        "kspacing": 0.14,
        "ecutwfc": 100.0,
        "scf_thr": 1e-7,
        "force_thr_ev": 0.01,
        "calculation": "scf",
        "basis_type": "lcao",
        "dft_functional": "pbe",
        "mixing_beta": 0.8,
        "smearing_sigma": 0.015,
        "smearing_method": "gauss",
        "relax_method": "cg",
        "relax_nmax": 60,
        "ks_solver": "genelpa",
    },
    "plotting": {
        "style": "abacuscopilot",
        "dpi": 300,
        "figure_format": "png",
        "figure_size": [8, 6],
        "color_cycle": "tab10",
        "font_size": 12,
        "show_fermi": True,
    },
    "paths": {
        "abacus_binary": "abacus",
        "mpirun": "mpirun",
        "abacus_source": "",
        "slurm_env_file": "",
        "sub_script": "",
        "sub_script_dp": "",
        "deepmd_python": "",
        "abacus_dp_binary": "",
        "abacus_dp_version": "",
    },
    "libraries": {
        "pseudo_library": _detect_library_dirs(".upf"),
        "orbital_library": _detect_library_dirs(".orb"),
        "family": "sg15",
        "rcut_policy": "7",
        "families": {},
    },
}


def _valid_library_dirs(value: Any, ext: str) -> list[str]:
    """Return the subset of *value* (str or list) that are real library dirs.

    A dir is valid if it exists and contains at least one *{ext} file
    (recursively).  Resolves ``~`` and symlinks.
    """
    dirs = [value] if isinstance(value, str) and value else list(value or [])
    valid: list[str] = []
    for d in dirs:
        p = Path(d).expanduser()
        if p.is_dir() and any(
            f.is_file() and f.name.lower().endswith(ext) for f in p.rglob("*")
        ):
            valid.append(str(p.resolve()))
    return valid


def _get_config_dir() -> Path:
    """Get the abacuscopilot config directory (~/.abacuscopilot/)."""
    return Path.home() / ".abacuscopilot"


def _get_config_path() -> Path:
    """Get the full path to the config file."""
    return _get_config_dir() / "config.yaml"


def _ensure_config_dir() -> None:
    """Create the config directory if it doesn't exist."""
    _get_config_dir().mkdir(parents=True, exist_ok=True)




def load_config() -> dict[str, Any]:
    """Load configuration from ~/.abacuscopilot/config.yaml.

    If the file doesn't exist, tries to migrate from the legacy
    ~/.abacuskit/config.yaml location first.  Falls back to a fresh
    default config when neither exists.

    Returns:
        Configuration dictionary.
    """
    config_path = _get_config_path()
    legacy_path = Path.home() / ".abacuskit" / "config.yaml"

    if not config_path.exists():
        _ensure_config_dir()
        if legacy_path.exists():
            # Migrate legacy config, updating paths in case the project dir was renamed
            with open(legacy_path) as f:
                user_config = yaml.safe_load(f) or {}
            save_config(user_config)
            return _deep_merge(DEFAULT_CONFIG, user_config)
        save_config(DEFAULT_CONFIG)
        return dict(DEFAULT_CONFIG)

    with open(config_path) as f:
        user_config = yaml.safe_load(f) or {}

    # Deep merge with defaults for any missing keys
    merged = _deep_merge(DEFAULT_CONFIG, user_config)

    # Sanitize library paths.  A configured list is kept verbatim for as long as
    # every entry still exists; the moment one is gone -- a library folder
    # moved, as in the PP-Orb/ reorg that nested each series under one folder --
    # the bundled roots are re-detected and merged in.  The old rule (keep the
    # list whenever *any* entry survived) silently lost a whole series: a config
    # listing the old flat SG15 path beside the lanthanide path kept only the
    # latter, and then every SG15 element failed to resolve.
    #
    # The re-detected roots are merged with the survivors rather than written
    # over them, because replacing the list would drop a `family: custom` user's
    # own hand-listed directories -- the same silent-wrong-pseudopotential
    # failure this whole block guards against.
    #
    # Two exclusions, doing two different jobs: _external_series_tops() matches
    # by folder NAME and keeps a known series such as Dojo-NC-FR out even on a
    # machine that never registered it; _registered_family_dirs() matches by
    # PATH and keeps a registered family's dirs out wherever they happen to live.
    libs = merged.setdefault("libraries", {})
    exclude = _registered_family_dirs(merged)
    for key, ext in (("pseudo_library", ".upf"), ("orbital_library", ".orb")):
        value = libs.get(key, "")
        configured = [value] if isinstance(value, str) and value else list(value or [])
        valid = _valid_library_dirs(value, ext)
        if valid and len(valid) == len(configured):
            libs[key] = valid          # nothing configured is gone — leave it be
            continue
        detected = _detect_library_dirs(ext, exclude_containing=exclude)
        libs[key] = detected + [d for d in valid if d not in detected]

    # The active family owns the library lists, and the default family is SG15.
    # A persisted list can still hold entries from a series picked earlier (or
    # hand-edited in), and the resolver searches every listed dir globally, so
    # `family: sg15` beside a Dojo entry resolved a Si job to Dojo `Si.upf` +
    # its SZ `..._1s1p.orb` instead of SG15's `Si_ONCV_PBE-1.0.upf` +
    # `2s2p1d`.  Re-derive the bundled SG15 (+ lanthanide) roots whenever SG15
    # is the active family; an external family's own dirs are written by
    # materialize_family() when it is picked, and `custom` keeps its
    # hand-configured lists.  Replace only when detection finds something, so
    # a machine without PP-Orb/ keeps whatever it had.
    if (libs.get("family") or "sg15") == "sg15":
        for key, ext in (("pseudo_library", ".upf"), ("orbital_library", ".orb")):
            detected = _detect_library_dirs(ext, exclude_containing=exclude)
            if detected:
                libs[key] = detected

    # A `defaults.pseudo_dir` / `orbital_dir` pointing at a directory that no
    # longer exists is a leftover from an older wizard run; ABACUS would fail
    # there.  "./" is the value that means "the job directory, where the
    # matching files are copied" — which is what the auto-copy step produces.
    for key in ("pseudo_dir", "orbital_dir"):
        _reset_missing_default_dir(merged, key)

    return merged


def _reset_missing_default_dir(merged: dict, key: str) -> bool:
    """Send ``defaults.<key>`` back to ``"./"`` when its directory is gone.

    Returns True if the value was changed.  A blank value, ``"./"`` itself, and
    a path that still exists are all left alone.
    """
    defaults = merged.setdefault("defaults", {})
    value = defaults.get(key)
    if not isinstance(value, str) or not value or value in ("./", "."):
        return False
    if Path(value).expanduser().exists():
        return False
    defaults[key] = "./"
    return True


#: One-line note written to the RIGHT of a value, in the style of VASPKIT's
#: ``~/.vaspkit``.  Keep every entry to a single short line.
#:
#: The layout is the point, so it is worth saying what it is not.  A first
#: attempt put multi-line explanations ABOVE each key: the file roughly doubled
#: in length, and since comments were then the majority of the lines it read as
#: a wall of prose with the values -- the only thing anyone opens the file for --
#: buried in it.  A note beside the value costs no extra line and stays out of
#: the way.  Anything that needs a paragraph goes in USER_GUIDE.md section 20.
_FIELD_HELP: dict[str, str] = {
    # defaults -- plain ABACUS keywords
    "defaults.pseudo_dir": "pseudopotential dir; './' = job dir (matching files are copied in)",
    "defaults.orbital_dir": "orbital dir; './' = job dir (matching files are copied in)",
    "defaults.kspacing": "k-spacing in 1/Bohr; ignored when a KPT file is present",
    "defaults.ecutwfc": "plane-wave cutoff, Ry",
    "defaults.scf_thr": "SCF convergence threshold",
    "defaults.force_thr_ev": "force convergence, eV/Angstrom",
    "defaults.calculation": "default task: scf | relax | cell-relax | md | ...",
    "defaults.basis_type": "lcao (numeric orbitals) | pw (plane waves)",
    "defaults.dft_functional": "exchange-correlation functional",
    "defaults.mixing_beta": "charge mixing; lower it (0.2) when SCF will not converge",
    "defaults.smearing_sigma": "smearing width, Ry",
    "defaults.smearing_method": "smearing scheme",
    "defaults.relax_method": "relaxation algorithm",
    "defaults.relax_nmax": "maximum ionic steps",
    "defaults.ks_solver": "eigenvalue solver",
    # plotting
    "plotting.style": "figure style; names are defined in plotting/style.py",
    "plotting.dpi": "figure resolution",
    "plotting.figure_format": "png | pdf | eps | svg | ...",
    "plotting.figure_size": "figure size in inches, [width, height]",
    "plotting.color_cycle": "NOT USED YET -- nothing in the code reads this",
    "plotting.font_size": "base font size",
    "plotting.show_fermi": "NOT USED YET -- nothing in the code reads this",
    # paths
    "paths.abacus_binary": "ABACUS executable; a bare name is looked up on $PATH",
    "paths.mpirun": "MPI launcher, used for parallel runs",
    "paths.abacus_source": "ABACUS source tree, for the abacuslite Python interface",
    "paths.slurm_env_file": "sourced at the top of every SLURM job",
    "paths.sub_script": "sbatch template, copied into the job dir",
    "paths.sub_script_dp": "as above for ABACUS-DP jobs; empty or missing = use sub_script",
    "paths.deepmd_python": "NOT USED YET -- nothing in the code reads this",
    "paths.abacus_dp_binary": "ABACUS-DP binary; empty = probe abacus_binary",
    "paths.abacus_dp_version": "DeepMD-kit version it was built against -- NOT ABACUS's",
    # libraries
    "libraries.pseudo_library": "searched in order; a stale entry re-detects from PP-Orb/",
    "libraries.orbital_library": "as above, for *.orb",
    "libraries.family": "sg15 | apns | apns/<variant> | dojoncfr[/<tier>] | custom",
    "libraries.rcut_policy": "'7' | 'min' | 'max'; dojoncfr only -- quote it",
    "libraries.families": "external series you downloaded yourself, keyed by id",
}

#: One line above a section heading.  VASPKIT groups its keys the same way.
_SECTION_HELP: dict[str, str] = {
    "defaults": "Defaults for generated ABACUS INPUT files -- plain ABACUS keywords.",
    "plotting": "Figure defaults.",
    "paths": "External programs abacuscopilot runs.",
    "libraries": "Pseudopotentials and orbitals.",
}

_CONFIG_HEADER = (
    "# abacuscopilot configuration\n"
    "#\n"
    "# Plain YAML -- edit any value in place. The notes on the right are for\n"
    "# reference only; deleting them changes nothing. Longer descriptions and\n"
    "# the meaning of every key: USER_GUIDE.md, section 20.\n"
    "#\n"
    "# Generated files are only written when the file is missing or you run the\n"
    "# wizard (task 9901); reading the config never rewrites it.\n"
)

#: Indent of a key inside its section.
_INDENT = "  "

#: How far right a note may push the value before it is put on the next line
#: instead -- past this the eye can no longer connect the two.
_MAX_GAP = 62

#: Widest a value may be and still be written inline (``[8, 6]`` rather than a
#: two-line block sequence).
_MAX_INLINE = 40


def _yaml(obj: Any, flow: bool) -> str:
    """Dump *obj* without the trailing newline or ``...`` end marker."""
    text = yaml.safe_dump(
        obj,
        default_flow_style=flow,
        sort_keys=False,
        allow_unicode=True,
        # Don't fold long scalars across lines: a path is the long value here,
        # and a wrapped one is hard to both read and hand-edit.
        width=10_000,
    ).rstrip("\n")
    return text[:-4] if text.endswith("\n...") else text


def _dump_pair(key: str, value: Any) -> list[str]:
    """Render one ``key: value`` pair as lines, unpadded and unindented."""
    block = _yaml({key: value}, flow=False).split("\n")
    # A short collection reads better inline, where it can share its line with
    # the note -- ``figure_size: [8, 6]``.  A long one must not: pseudo_library
    # holds full paths, and inline they would run off the side of the terminal.
    if len(block) > 1 and isinstance(value, (list, dict)) and value:
        inline = _yaml(value, flow=True)
        if "\n" not in inline and len(f"{key}: {inline}") <= _MAX_INLINE:
            return [f"{key}: {inline}"]
    return block


def render_config(config: dict[str, Any]) -> str:
    """Serialize *config* with each value's note aligned to its right.

    Values are dumped by PyYAML, one pair at a time, so ``libraries.families``
    (an arbitrary user-defined tree) round-trips exactly as it would from a
    plain ``safe_dump``; only the notes and their column are ours.  A value
    spanning several lines -- a list, a nested mapping -- takes its note on the
    key's line, and its continuation lines stay uncommented.
    """
    blocks: list[tuple[str, list[tuple[str, list[str]]]]] = [
        (section, [(key, _dump_pair(key, value)) for key, value in values.items()])
        for section, values in config.items() if isinstance(values, dict) and values
    ]
    # One column for the whole file, not one per section: a single straight
    # margin is what the eye follows, and per-section columns put four of them
    # in a forty-line file.  Capped, because past the cap the note is too far
    # from its value to be read as belonging to it.
    column = min(
        max((len(_INDENT) + len(lines[0]) for _, pairs in blocks for _, lines in pairs),
            default=0) + 2,
        len(_INDENT) + _MAX_GAP,
    )

    out: list[str] = [_CONFIG_HEADER]
    for section, values in config.items():
        out.append("\n")
        heading = _SECTION_HELP.get(section)
        if heading:
            out.append(f"# {heading}\n")
        if not isinstance(values, dict):
            out.append("\n".join(_dump_pair(section, values)) + "\n")
            continue
        if not values:
            # A bare ``paths:`` loads back as None, not {}.
            out.append(f"{section}: {{}}\n")
            continue

        out.append(f"{section}:\n")
        for key, lines in dict(blocks)[section]:
            head = _INDENT + lines[0]
            note = _FIELD_HELP.get(f"{section}.{key}")
            if note:
                head += " " * max(2, column - len(head)) + f"# {note}"
            out.append(head + "\n")
            for cont in lines[1:]:
                # A block sequence is the one continuation PyYAML does not
                # indent, leaving "- item" in the key's own column where it
                # reads as a sibling.  A nested mapping already arrives
                # indented, and indenting it again would over-nest it.
                extra = "  " if cont.startswith("- ") else ""
                out.append(_INDENT + extra + cont + "\n")
    return "".join(out)


def save_config(config: dict[str, Any]) -> None:
    """Save configuration to ~/.abacuscopilot/config.yaml.

    Args:
        config: Configuration dictionary to save.
    """
    _ensure_config_dir()
    config_path = _get_config_path()
    with open(config_path, "w") as f:
        f.write(render_config(config))


def get_config_value(key_path: str, default=None) -> Any:
    """Get a specific config value using dot-separated path.

    Example:
        get_config_value("defaults.kspacing") -> 0.04

    Args:
        key_path: Dot-separated key path (e.g., "defaults.kspacing").
        default: Value to return if key not found.

    Returns:
        The config value, or default if not found.
    """
    config = load_config()
    keys = key_path.split(".")
    value = config
    for key in keys:
        if isinstance(value, dict):
            value = value.get(key)
            if value is None:
                return default
        else:
            return default
    return value


def _deep_merge(base: dict, override: dict) -> dict:
    """Deep merge two dictionaries, with override taking precedence.

    Nested dicts are merged recursively. Non-dict values from override
    replace those in base.
    """
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result
