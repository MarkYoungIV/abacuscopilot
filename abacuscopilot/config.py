"""Configuration management for abacuscopilot.

Reads and writes YAML config from ~/.abacuscopilot/config.yaml.
Settings include default paths, pseudopotential directories,
plotting preferences, and user-defined presets.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def _detect_library_dirs(ext: str) -> list[str]:
    """Auto-detect library directories inside the package root.

    Looks for top-level library folders under ``PP-Orb/`` (preferred) or the
    legacy ``Pseudopotential/`` / ``Orbitals/`` directories, and returns those
    roots that contain at least one *{ext} file (searched recursively).

    A folder is a library root if it directly or indirectly holds the actual
    files, e.g.:
      PP-Orb/SG15-Version1p0_Pseudopotential/           (has .upf below)
      PP-Orb/SG15-Version1p0__StandardOrbitals-Version2p0/  (has .orb)
      PP-Orb/lanthanides-f--core.icmod1/                (has .UPF *and* .orb)

    Returns a list of absolute paths (possibly empty).
    """
    pkg_root = Path(__file__).resolve().parent.parent
    roots: list[Path] = []
    pp_orb = pkg_root / "PP-Orb"
    if pp_orb.is_dir():
        roots = [d for d in sorted(pp_orb.iterdir())
                 if d.is_dir() and not d.name.startswith(".")]
    else:
        for name in ("Pseudopotential", "Orbitals"):
            d = pkg_root / name
            if d.is_dir():
                roots.append(d)

    found: list[str] = []
    for root in roots:
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
        "abacus_source": "",     # ABACUS source tree (for abacuslite PYTHONPATH)
        "slurm_env_file": "",    # shell script to source in SLURM jobs (CUDA, compiler, etc.)
        "sub_script": "",        # Slurm sbatch template — copied alongside INPUT (101-110)
        "sub_script_dp": "",     # Slurm sbatch template for ABACUS-DP (113)
        "deepmd_python": "",     # Python binary with deepmd-kit (for Deep Potential batch force calculation)
    },
    "libraries": {
        "pseudo_library": _detect_library_dirs(".upf"),
        "orbital_library": _detect_library_dirs(".orb"),
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

    # Sanitize library paths: drop stale/empty dirs (e.g. the old
    # Pseudopotential/Orbitals locations after a reorg), fall back to the
    # auto-detected defaults.  Result values are always lists of valid dirs.
    libs = merged.setdefault("libraries", {})
    for key, ext in (("pseudo_library", ".upf"), ("orbital_library", ".orb")):
        valid = _valid_library_dirs(libs.get(key, ""), ext)
        libs[key] = valid if valid else DEFAULT_CONFIG["libraries"].get(key, [])
    return merged


def save_config(config: dict[str, Any]) -> None:
    """Save configuration to ~/.abacuscopilot/config.yaml.

    Args:
        config: Configuration dictionary to save.
    """
    _ensure_config_dir()
    config_path = _get_config_path()
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, default_flow_style=False, sort_keys=False)


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
