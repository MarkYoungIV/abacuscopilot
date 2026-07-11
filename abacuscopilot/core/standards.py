"""Standard output规范 — single authoritative source for basis_type / ks_solver.

This module encodes AbacusCopilot's standard output规范, the rules that every
task must follow so generated INPUT/STRU files are internally consistent:

- basis_type == "pw":
    * ks_solver = dav_subspace (on CPU *and* GPU — PW has no GPU eigensolver split)
    * STRU writes ATOMIC_SPECIES with .upf only; NO NUMERICAL_ORBITAL section
    * only .upf pseudopotentials are copied
- basis_type == "lcao" (or "lcao_in_pw"):
    * ks_solver = genelpa on CPU, cusolver on a single GPU card
    * STRU writes both .upf and NUMERICAL_ORBITAL (.orb) sections
    * both .upf and .orb files are copied

Authoritative reference: ABACUS user manual (device/ks_solver section) and the
MD tutorial note "ks_solver 改为 genelpa, gpu 单卡用 cusolver".

Keep this module dependency-free (no config / task imports) to avoid import
cycles — it is pure logic that any layer may call.
"""

from __future__ import annotations

# Valid eigensolvers per basis type (for validating user-supplied --solver).
_ALLOWED_SOLVERS = {
    "pw": {"dav_subspace", "cg", "bpcg", "dav"},
    "lcao": {"genelpa", "cusolver", "cusolvermp", "elpa", "scalapack_gvx", "lapack"},
}


def is_lcao_basis(basis_type: str | None) -> bool:
    """True if the basis type is any LCAO variant (lcao, lcao_in_pw)."""
    return bool(basis_type) and str(basis_type).startswith("lcao")


def solver_for(basis_type: str | None, device: str = "cpu") -> str:
    """Authoritative (basis_type, device) -> ks_solver mapping.

    pw            -> dav_subspace   (regardless of device)
    lcao on cpu   -> genelpa
    lcao on gpu   -> cusolver       (single GPU card)
    """
    if not is_lcao_basis(basis_type):
        return "dav_subspace"
    return "cusolver" if str(device).lower() == "gpu" else "genelpa"


def validate_solver(
    basis_type: str | None,
    ks_solver: str | None,
    device: str = "cpu",
) -> tuple[str, str | None]:
    """Validate a user-supplied ks_solver against the standard规范.

    Returns (corrected_solver, warning_or_None):
    - If ks_solver is empty/None, returns the规范 default for (basis_type, device).
    - If ks_solver is valid for the basis type, returns it unchanged.
    - If ks_solver is invalid for the basis type, returns the规范 default plus a
      warning message describing the auto-correction.
    """
    expected = solver_for(basis_type, device)
    if not ks_solver:
        return expected, None
    key = "lcao" if is_lcao_basis(basis_type) else "pw"
    if ks_solver not in _ALLOWED_SOLVERS[key]:
        return (
            expected,
            f"ks_solver '{ks_solver}' is invalid for basis_type "
            f"'{basis_type}'; corrected to '{expected}'.",
        )
    return ks_solver, None
