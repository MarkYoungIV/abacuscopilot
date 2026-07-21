"""Configuration management for abacuscopilot.

Reads and writes YAML config from ~/.abacuscopilot/config.yaml.
Settings include default paths, pseudopotential directories,
plotting preferences, and user-defined presets.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def _auto_detect_lib(subdir_name: str) -> str:
    """Auto-detect a library directory inside the package root.

    Finds the subdirectory containing actual .upf / .orb files
    (recursively), or returns '' if none found.
    """
    try:
        pkg_root = Path(__file__).resolve().parent.parent
        lib_dir = pkg_root / subdir_name
        if lib_dir.is_dir():
            ext = ".upf" if subdir_name == "Pseudopotential" else ".orb"
            for d in sorted(lib_dir.rglob("*"), reverse=True):
                if d.is_dir() and not d.name.startswith("."):
                    if list(d.glob(f"*{ext}")):
                        return str(d.resolve())
    except Exception:
        pass
    return ""


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
        "deepmd_python": "",     # Python binary with deepmd-kit (for 1511 batch force calc)
    },
    "libraries": {
        "pseudo_library": _auto_detect_lib("Pseudopotential"),
        "orbital_library": _auto_detect_lib("Orbitals"),
    },
}


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
    return _deep_merge(DEFAULT_CONFIG, user_config)


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
