"""Reader and writer for ABACUS INPUT files.

The INPUT file is the master control file for ABACUS calculations.
It contains keyword-value pairs specifying all calculation parameters.

Format:
    INPUT_PARAMETERS
    keyword1    value1
    keyword2    value2
    ...

Lines starting with # or / are comments. Parameters can be in any order.
Boolean values accept: True/False, 1/0, T/F (case-insensitive).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from abacuscopilot.core.exceptions import FileFormatError, FileNotFoundError_
from abacuscopilot.core.models import InputParams

# Regex for parsing keyword-value lines
# Matches: keyword = value, or keyword value
_KV_RE = re.compile(
    r"^\s*([a-zA-Z_][a-zA-Z0-9_]*)\s+(.+?)\s*(?:#.*)?$"
)

# Valid boolean representations
_TRUE_VALUES = {"true", "t", "1", "yes", "y"}
_FALSE_VALUES = {"false", "f", "0", "no", "n"}


def _parse_value(raw: str) -> Any:
    """Parse a raw string value into the appropriate Python type.

    Tries in order: bool, int, float, string.
    """
    raw = raw.strip()

    # Boolean
    low = raw.lower()
    if low in _TRUE_VALUES:
        return True
    if low in _FALSE_VALUES:
        return False

    # Integer
    try:
        return int(raw)
    except ValueError:
        pass

    # Float
    try:
        return float(raw)
    except ValueError:
        pass

    # String (strip quotes if present)
    if (raw.startswith('"') and raw.endswith('"')) or \
       (raw.startswith("'") and raw.endswith("'")):
        return raw[1:-1]

    return raw


def _format_value(value: Any) -> str:
    """Format a Python value for writing to an INPUT file.

    Booleans -> 1/0, floats -> .8g, strings -> as-is (quoted if needed).
    """
    if isinstance(value, bool):
        return "1" if value else "0"
    elif isinstance(value, (int, float)):
        if isinstance(value, float):
            return f"{value:.8g}"
        return str(value)
    elif isinstance(value, list):
        return " ".join(str(v) for v in value)
    else:
        return str(value)


def read_input(filepath: str | Path) -> InputParams:
    """Read and parse an ABACUS INPUT file.

    Args:
        filepath: Path to the INPUT file.

    Returns:
        InputParams object with all parsed parameters.

    Raises:
        FileNotFoundError_: If the file doesn't exist.
        FileFormatError: If the file doesn't start with INPUT_PARAMETERS.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError_(str(filepath), "INPUT file not found")

    with open(filepath) as f:
        lines = f.readlines()

    if not lines:
        raise FileFormatError(str(filepath), "File is empty")

    # Check for INPUT_PARAMETERS header
    first_line = lines[0].strip()
    if not first_line.startswith("INPUT_PARAMETERS"):
        raise FileFormatError(
            str(filepath),
            f"Expected 'INPUT_PARAMETERS' on first line, got '{first_line}'"
        )

    params = InputParams()
    for lineno, line in enumerate(lines[1:], start=2):
        stripped = line.strip()

        # Skip empty lines and comments
        if not stripped or stripped.startswith("#") or stripped.startswith("/"):
            continue

        # Parse keyword = value
        match = _KV_RE.match(line)
        if not match:
            continue  # Skip unparseable lines

        keyword = match.group(1).lower()
        raw_value = match.group(2).strip()

        # Remove inline comment from value
        if "#" in raw_value:
            raw_value = raw_value.split("#")[0].strip()

        value = _parse_value(raw_value)

        # Route to known field or extras
        if keyword in params._known_params:
            setattr(params, keyword, value)
        else:
            params.extras[keyword] = value

    return params


def write_input(
    params: InputParams,
    filepath: str | Path = "INPUT",
    comment: str = "",
) -> None:
    """Write an InputParams object to an ABACUS INPUT file.

    Parameters are written in groups with section comments for readability.

    Args:
        params: InputParams object to write.
        filepath: Output file path.
        comment: Optional comment to add after INPUT_PARAMETERS line.
    """
    filepath = Path(filepath)

    # Map Python field names → ABACUS INPUT keyword names (for names that differ)
    # Most fields match the ABACUS keyword directly, but a few need explicit mapping.
    _keyword_map: dict[str, str] = {
        # force_thr_ev is the canonical ABACUS keyword; force_thr is a model alias
        "force_thr": "force_thr_ev",
    }

    # Define parameter groups for organized output.
    groups: list[tuple[str, list[str]]] = [
        ("System variables", [
            "calculation", "symmetry", "kspacing", "gamma_only", "precision",
        ]),
        ("Plane wave related variables", [
            "ecutwfc", "ecutrho", "pw_diag_nmax", "pw_diag_ndim",
        ]),
        ("Electronic structure", [
            "basis_type", "ks_solver", "smearing_method", "smearing_sigma",
            "mixing_type", "mixing_beta", "scf_nmax", "scf_thr",
            "nbands", "nelec", "dft_functional", "nspin", "noncolin",
            "lspinorb",
        ]),
        ("SCF control", [
            "chg_extrap",
        ]),
        ("Geometry relaxation", [
            "relax_method", "relax_nmax",
            "cal_force", "force_thr_ev", "cal_stress", "stress_thr",
            "fixed_axes",
        ]),
        ("I/O", [
            "stru_file", "kpoint_file", "pseudo_dir", "orbital_dir",
            "read_file_dir", "restart_load",
        ]),
        ("MD ensemble control", [
            "md_type", "md_pmode", "press1", "press2", "press3",
        ]),
        ("MD process control", [
            "md_nstep", "md_dt", "md_tfirst", "md_tlast", "md_damp",
        ]),
        ("DP model", [
            "esolver_type", "pot_file",
        ]),
        ("Output set", [
            "md_restartfreq", "md_dumpfreq", "out_level",
            "dump_force", "dump_vel",
        ]),
        ("Output", [
            "out_chg", "out_pot", "out_dos", "out_band",
            "out_proj_band", "out_bandgap", "out_wfc_pw",
            "out_wfc_lcao", "out_wfc_r", "out_elf",
            "out_mat_hs", "out_mat_hs2", "out_dm", "out_mul",
        ]),
        ("DOS", [
            "dos_emin_ev", "dos_emax_ev", "dos_edelta_ev", "dos_sigma",
        ]),
    ]

    # Core physics groups: write every parameter the template explicitly sets.
    _CORE_GROUPS = {"System variables", "Plane wave related variables",
                    "Electronic structure", "SCF control", "Geometry relaxation",
                    "MD process control"}

    def _is_explicitly_set(val: Any) -> bool:
        """Return True if *val* looks like an explicitly-set parameter."""
        if val is None:
            return False
        if isinstance(val, str) and val == "":
            return False
        return True

    # Inline comments: maps param_name → comment text (without leading "# ")
    _inline_comments: dict[str, str] = params.extras.get("_inline_comments", {})

    # Section-level commented hints (e.g. #kspacing in System variables)
    _section_hints: dict[str, list[str]] = params.extras.get("_section_hints", {})

    # Bottom-of-file commented hints (e.g. #cal_force, #cal_stress)
    _comment_hints: dict[str, str] = params.extras.get("_comment_hints", {})

    with open(filepath, "w") as f:
        # Header
        if comment:
            f.write(f"INPUT_PARAMETERS  # {comment}\n")
        else:
            f.write("INPUT_PARAMETERS\n")

        # Dedup: if force_thr_ev is at default, use force_thr value
        if getattr(params, "force_thr_ev", None) == getattr(InputParams(), "force_thr_ev", None):
            params.force_thr_ev = params.force_thr

        written: set[str] = set()
        for group_name, keys in groups:
            group_params = {}
            for key in keys:
                if key == "force_thr" and "force_thr_ev" in keys:
                    continue
                val = getattr(params, key, None)
                if val is None:
                    continue

                abacus_key = _keyword_map.get(key, key)

                if group_name in _CORE_GROUPS:
                    # Core physics groups: write if value was explicitly set
                    # AND appears in the template's touched-keys (or no template used)
                    if _is_explicitly_set(val):
                        template_keys = params.extras.get("_template_keys", None)
                        if template_keys is None or key in template_keys:
                            group_params[abacus_key] = val
                            written.add(key)
                else:
                    # Conditional groups: write when value differs from default,
                    # OR when explicitly touched by a template
                    default_val = getattr(InputParams(), key, None)
                    template_keys = params.extras.get("_template_keys", None)
                    if val != default_val or (template_keys is not None and key in template_keys):
                        group_params[abacus_key] = val
                        written.add(key)

            if group_params or group_name in _section_hints:
                f.write(f"\n# {group_name}\n")
                # Section-level commented hints first (e.g. #kspacing 0.14)
                for hint_line in _section_hints.get(group_name, []):
                    f.write(f"#{hint_line}\n")
                for abacus_key, val in group_params.items():
                    val_str = _format_value(val)
                    line = f"{abacus_key:<20s} {val_str}"
                    if abacus_key in _inline_comments:
                        line += f"  # {_inline_comments[abacus_key]}"
                    f.write(line + "\n")

        # Write any remaining known params not in groups (conditional I/O-style check)
        remaining_known = []
        for key in params._known_params:
            if key not in written and key not in ("force_thr", "_known_params"):
                val = getattr(params, key, None)
                if val is not None and _is_explicitly_set(val):
                    default_val = getattr(InputParams(), key, None)
                    if val != default_val:
                        abacus_key = _keyword_map.get(key, key)
                        remaining_known.append((abacus_key, val))
                        written.add(key)

        if remaining_known:
            f.write("\n# Other Parameters\n")
            for abacus_key, val in remaining_known:
                val_str = _format_value(val)
                line = f"{abacus_key:<20s} {val_str}"
                if abacus_key in _inline_comments:
                    line += f"  # {_inline_comments[abacus_key]}"
                f.write(line + "\n")

        # Write bottom-of-file commented hints (e.g. #cal_force, #cal_stress)
        if _comment_hints:
            f.write("\n# Other Parameters\n")
            for key, hint_text in _comment_hints.items():
                f.write(f"#{key:<19s} {hint_text}\n")

        # Write extras (skip internal meta-keys)
        _meta_keys = {"_template_keys", "_comment_hints", "_inline_comments", "_section_hints"}
        user_extras = {k: v for k, v in params.extras.items() if k not in _meta_keys}
        if user_extras:
            f.write("\n# Additional\n")
            for key, val in user_extras.items():
                f.write(f"{key:<20s} {_format_value(val)}\n")


def validate_input(params: InputParams) -> list[str]:
    """Validate InputParams for common issues.

    Args:
        params: InputParams object to validate.

    Returns:
        List of warning/error messages. Empty list means no issues found.
    """
    warnings = []

    # Check required params
    if params.ntype <= 0:
        warnings.append("ntype must be > 0")

    if params.calculation not in (
        "scf", "nscf", "relax", "cell-relax", "md",
        "get_pchg", "get_wf", "get_s", "gen_bessel",
        "gen_opt_abfs", "test_memory", "test_neighbour",
    ):
        warnings.append(f"Unknown calculation type: {params.calculation}")

    if params.basis_type not in ("pw", "lcao", "lcao_in_pw"):
        warnings.append(f"Unknown basis_type: {params.basis_type}")

    if params.ecutwfc <= 0:
        warnings.append("ecutwfc must be > 0")

    if params.scf_thr <= 0:
        warnings.append("scf_thr must be > 0")

    # LCAO-specific checks
    if params.basis_type == "lcao":
        if params.ks_solver not in (
            "genelpa", "scalapack_gvx", "lapack", "cusolver", "elpa",
        ):
            warnings.append(f"Unknown LCAO ks_solver: {params.ks_solver}")

    return warnings
