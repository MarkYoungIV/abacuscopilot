"""Core data models for abacuscopilot.

These dataclasses are the canonical in-memory representation of ABACUS
structures, k-points, and input parameters used throughout the codebase.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from ase import Atoms


# =============================================================================
# Lattice
# =============================================================================


@dataclass
class Lattice:
    """Crystal lattice with vectors and scaling constant.

    The actual lattice vectors in Bohr are: cell = constant * vectors.

    Attributes:
        constant: Scaling factor in Bohr (the LATTICE_CONSTANT value).
        vectors: 3x3 matrix of unitless lattice vectors (rows).
    """

    constant: float = 1.0
    vectors: np.ndarray = field(
        default_factory=lambda: np.eye(3)
    )

    def __post_init__(self):
        if self.vectors.shape != (3, 3):
            raise ValueError(f"Lattice vectors must be 3x3, got {self.vectors.shape}")

    @property
    def cell(self) -> np.ndarray:
        """Actual cell vectors in Bohr: constant * vectors."""
        return self.constant * self.vectors

    @property
    def cell_angstrom(self) -> np.ndarray:
        """Actual cell vectors in Angstrom."""
        from abacuscopilot.core.constants import BOHR_TO_ANGSTROM
        return self.cell * BOHR_TO_ANGSTROM

    @property
    def volume(self) -> float:
        """Cell volume in Bohr^3."""
        return abs(np.linalg.det(self.cell))

    @property
    def volume_angstrom(self) -> float:
        """Cell volume in Angstrom^3."""
        from abacuscopilot.core.constants import BOHR_TO_ANGSTROM
        return self.volume * (BOHR_TO_ANGSTROM ** 3)

    @property
    def reciprocal_cell(self) -> np.ndarray:
        """Reciprocal lattice vectors in 2π/Bohr (rows = b1, b2, b3)."""
        return 2.0 * np.pi * np.linalg.inv(self.cell).T

    @classmethod
    def from_cell(cls, cell: np.ndarray, constant: float = 1.0) -> "Lattice":
        """Create Lattice from actual cell vectors (Bohr)."""
        return cls(constant=constant, vectors=cell / constant)

    @classmethod
    def from_cell_angstrom(cls, cell: np.ndarray, constant: float | None = None) -> "Lattice":
        """Create Lattice from cell vectors in Angstrom.

        By default, uses the ABACUS convention: LATTICE_CONSTANT = 1.889726 (Bohr)
        and vectors equal to the Angstrom values. This keeps vector numbers readable
        (same as in the CIF/POSCAR) while producing the correct cell in Bohr:
            actual cell = constant * vectors = 1.889726 * cell_ang → Bohr.
        """
        from abacuscopilot.core.constants import ANGSTROM_TO_BOHR
        if constant is None:
            constant = ANGSTROM_TO_BOHR
        return cls(constant=constant, vectors=np.copy(cell))


# =============================================================================
# Atom
# =============================================================================


@dataclass
class Atom:
    """A single atom in a crystal structure.

    Attributes:
        species: Element label (e.g., 'Si', 'O').
        position: 3D position vector in the coordinate system specified
                  by the parent Structure.
        fix: Whether each Cartesian direction is free (True) or fixed (False)
             during relaxation. Defaults to all free (True, True, True).
        magmom: Initial magnetic moment for this atom.
        velocity: Initial velocity vector for MD simulations.
        angle1: Non-collinear spin angle from z-axis (degrees).
        angle2: Non-collinear spin angle in ab-plane (degrees).
    """

    species: str
    position: np.ndarray  # shape (3,)
    fix: tuple[bool, bool, bool] = (True, True, True)
    magmom: float = 0.0
    velocity: np.ndarray | None = None
    angle1: float | None = None
    angle2: float | None = None

    def __post_init__(self):
        if len(self.position) != 3:
            raise ValueError(f"Position must be length 3, got {len(self.position)}")


# =============================================================================
# Structure
# =============================================================================


@dataclass
class Structure:
    """Complete crystal structure as read from or written to an ABACUS STRU file.

    Attributes:
        lattice: Lattice object with constant and vectors.
        atoms: List of Atom objects.
        pseudo_files: Mapping from species label to pseudopotential file path.
        orbital_files: Mapping from species label to numerical orbital file path.
        coordinate_type: 'Direct', 'Cartesian', 'Cartesian_au',
                         'Cartesian_angstrom', etc.
        species_order: Ordered list of unique species (matches ATOMIC_SPECIES).
        magnetism: Default magnetism per species (for initial magmoms).
    """

    lattice: Lattice = field(default_factory=Lattice)
    atoms: list[Atom] = field(default_factory=list)
    pseudo_files: dict[str, str] = field(default_factory=dict)
    orbital_files: dict[str, str] = field(default_factory=dict)
    coordinate_type: str = "Direct"
    species_order: list[str] = field(default_factory=list)
    magnetism: dict[str, float] = field(default_factory=dict)

    @property
    def species(self) -> list[str]:
        """Return species labels for all atoms in order."""
        return [a.species for a in self.atoms]

    @property
    def num_atoms(self) -> int:
        """Total number of atoms."""
        return len(self.atoms)

    @property
    def num_species(self) -> int:
        """Number of unique species."""
        return len(self.species_order)

    @property
    def positions(self) -> np.ndarray:
        """Nx3 array of all atomic positions."""
        return np.array([a.position for a in self.atoms])

    def count_species(self, label: str) -> int:
        """Count atoms of a given species."""
        return sum(1 for a in self.atoms if a.species == label)

    def get_atoms_by_species(self, label: str) -> list[int]:
        """Get indices of atoms of a given species."""
        return [i for i, a in enumerate(self.atoms) if a.species == label]

    def to_ase(self) -> "Atoms":
        """Convert to ASE Atoms object.

        Requires ASE to be installed. Raises ASEImportError otherwise.
        """
        try:
            from ase import Atoms
        except ImportError as e:
            from abacuscopilot.core.exceptions import ASEImportError
            raise ASEImportError("Structure.to_ase()") from e


        cell_ang = self.lattice.cell_angstrom
        atoms = Atoms(
            symbols=self.species,
            positions=self.positions if self.coordinate_type in ("Cartesian", "Cartesian_angstrom")
                     else self.positions @ cell_ang,  # fractional -> Cartesian
            cell=cell_ang,
            pbc=True,
        )

        if self.coordinate_type.startswith("Direct"):
            atoms.set_scaled_positions(self.positions)

        return atoms

    @classmethod
    def from_ase(cls, atoms: "Atoms", **kwargs) -> "Structure":
        """Create Structure from ASE Atoms object.

        Requires ASE to be installed. Raises ASEImportError otherwise.
        """
        try:
            from ase import Atoms
        except ImportError as e:
            from abacuscopilot.core.exceptions import ASEImportError
            raise ASEImportError("Structure.from_ase()") from e


        lattice = Lattice.from_cell_angstrom(np.array(atoms.cell))
        structure = cls(
            lattice=lattice,
            coordinate_type="Direct",
            **kwargs,
        )

        # Determine species
        symbols = atoms.get_chemical_symbols()
        unique_species = []
        for s in symbols:
            if s not in unique_species:
                unique_species.append(s)
        structure.species_order = unique_species

        # Add atoms with fractional coordinates
        scaled_pos = atoms.get_scaled_positions()
        for i, (sym, pos) in enumerate(zip(symbols, scaled_pos)):
            structure.atoms.append(Atom(species=sym, position=pos))

        return structure


# =============================================================================
# KPoints
# =============================================================================


@dataclass
class KPoints:
    """K-point sampling specification.

    Supports all ABACUS KPT file modes:
    - Auto MP mesh (mode=0): grid and shift specified.
    - Explicit list (mode=Nk): list of k-point coordinates with weights.
    - Line mode: high-symmetry path for band structure.

    Attributes:
        mode: 'gamma', 'mp', 'direct', 'line', or 'line_cartesian'.
        grid: (N1, N2, N3) subdivisions for auto MP mesh.
        shift: (S1, S2, S3) shift for auto MP mesh.
        gamma_centered: Whether auto mesh is Gamma-centered (vs standard MP).
        explicit_kpoints: List of (kx, ky, kz, weight) for explicit mode.
        line_path: List of (start_xyz, end_xyz, npoints, label) for line mode.
                   The label is for the starting k-point; end label is from
                   the next segment's start.
        labels: Ordered list of high-symmetry point labels.
        label_positions: Cumulative k-distance of each label.
    """

    mode: str = "gamma"
    grid: tuple[int, int, int] | None = None
    shift: tuple[float, float, float] = (0.0, 0.0, 0.0)
    gamma_centered: bool = True
    explicit_kpoints: list[tuple[float, float, float, float]] = field(default_factory=list)
    line_path: list[dict] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)
    label_positions: list[float] = field(default_factory=list)

    def to_string(self) -> str:
        """Generate the KPT file content."""
        lines = ["K_POINTS"]

        if self.mode in ("gamma", "mp"):
            lines.append("0")
            lines.append("Gamma" if self.gamma_centered else "MP")
            g = self.grid or (1, 1, 1)
            s = self.shift
            lines.append(f"{g[0]} {g[1]} {g[2]}  {s[0]} {s[1]} {s[2]}")
            return "\n".join(lines) + "\n"

        elif self.mode == "direct":
            nk = len(self.explicit_kpoints)
            lines.append(str(nk))
            lines.append("Direct")
            for kx, ky, kz, w in self.explicit_kpoints:
                lines.append(f"  {kx:.10f}  {ky:.10f}  {kz:.10f}  {w:.10f}")
            return "\n".join(lines) + "\n"

        elif self.mode in ("line", "line_cartesian"):
            coord_type = "Line_Cartesian" if self.mode == "line_cartesian" else "Line"
            # Build the endpoint list from segments. A path may be discontinuous
            # (e.g. FCC: X->U then K->Gamma, where a segment's end != next
            # segment's start). At such a break the ABACUS convention (matching
            # abacustest) is to emit the break point with npoints=1, then start
            # the next branch fresh. Continuous joins share a single point.
            # Each entry: (xyz, npoints, label)
            points: list[tuple] = []
            segs = self.line_path
            if segs:
                first = segs[0]
                points.append((first["start"], first.get("npoints", 20),
                               first.get("label", "")))
                for idx, seg in enumerate(segs):
                    end = seg["end"]
                    end_label = seg.get("end_label", "")
                    is_last = idx == len(segs) - 1
                    if is_last:
                        # final endpoint terminates the path
                        points.append((end, 1, end_label))
                    else:
                        nxt = segs[idx + 1]
                        continuous = (tuple(end) == tuple(nxt["start"])
                                      and end_label == nxt.get("label", ""))
                        if continuous:
                            # shared point → carries next segment's npoints
                            points.append((end, nxt.get("npoints", 20), end_label))
                        else:
                            # discontinuity: close this branch (npoints=1),
                            # then open the next branch at its start
                            points.append((end, 1, end_label))
                            points.append((nxt["start"], nxt.get("npoints", 20),
                                           nxt.get("label", "")))
            lines.append(str(len(points)))
            lines.append(coord_type)
            for xyz, npts, label in points:
                label_str = f"  {label}" if label else ""
                lines.append(
                    f"  {xyz[0]:.10f}  {xyz[1]:.10f}  {xyz[2]:.10f}  {npts:>4d}{label_str}"
                )
            return "\n".join(lines) + "\n"

        raise ValueError(f"Unknown KPT mode: {self.mode}")

    def get_cumulative_distances(self, reciprocal_cell: np.ndarray) -> np.ndarray:
        """Compute cumulative k-distances along the line-mode path.

        Args:
            reciprocal_cell: 3x3 reciprocal lattice vectors (2π/Bohr, rows).

        Returns:
            Array of cumulative k-distances (len = total k-points).
        """
        if self.mode not in ("line", "line_cartesian"):
            raise ValueError("Cumulative distances only defined for line mode")

        segments = []
        for seg in self.line_path:
            start = np.array(seg["start"])
            end = np.array(seg["end"])
            npts = seg.get("npoints", 20)

            if self.mode == "line":
                # Fractional -> Cartesian in reciprocal space
                start_cart = start @ reciprocal_cell
                end_cart = end @ reciprocal_cell
            else:
                # line_cartesian: coordinates are already Cartesian
                start_cart = start
                end_cart = end

            seg_dist = np.linalg.norm(end_cart - start_cart)
            segments.append(seg_dist / npts * np.arange(npts))

        # Concatenate and compute cumulative sum
        segment_lengths = [np.linalg.norm(
            (np.array(s["end"]) - np.array(s["start"])) @ reciprocal_cell
            if self.mode == "line" else
            (np.array(s["end"]) - np.array(s["start"]))
        ) for s in self.line_path]

        # Build cumulative distances with labels at segment start points
        cumsum = 0.0
        all_dists = []
        self.label_positions = []
        self.labels = []

        # First label
        if self.line_path:
            first_seg = self.line_path[0]
            label = first_seg.get("label", "")
            self.labels.append(label)
            self.label_positions.append(0.0)

        for i, seg in enumerate(self.line_path):
            npts = seg.get("npoints", 20)
            seg_len = segment_lengths[i]
            dists = cumsum + np.linspace(0, seg_len, npts, endpoint=False)
            all_dists.append(dists)
            cumsum += seg_len

            # Label at end of this segment = start of next
            end_label = seg.get("end_label", "")
            if end_label:
                self.labels.append(end_label)
                self.label_positions.append(cumsum)

        return np.concatenate(all_dists) if all_dists else np.array([])


# =============================================================================
# InputParams
# =============================================================================


@dataclass
class InputParams:
    """All recognized parameters in an ABACUS INPUT file.

    Organized in logical groups matching the ABACUS documentation.
    Any unrecognized parameters are stored in `extras`.
    """

    # === System ===
    suffix: str = "ABACUS"
    ntype: int = 1
    calculation: str = "scf"
    esolver_type: str = "ksdft"
    pot_file: str = ""          # DP model path (only used when esolver_type = dp)
    symmetry: int = 1
    init_wfc: str = "atomic"
    init_chg: str = "atomic"
    kpar: int = 1
    bndpar: int = 1
    latname: str = "none"
    device: str = "cpu"
    kspacing: float = 0.0       # auto k-point spacing in 1/bohr (0 = disabled)
    precision: str = "double"    # single or double precision

    # === I/O ===
    stru_file: str = "STRU"
    kpoint_file: str = "KPT"
    pseudo_dir: str = "./"
    orbital_dir: str = "./"
    read_file_dir: str = ""
    restart_load: bool = False

    # === Plane-wave basis ===
    basis_type: str = "pw"
    ecutwfc: float = 100.0
    ecutrho: float = 0.0
    nx: int = 0
    ny: int = 0
    nz: int = 0
    pw_diag_thr: float = 0.01
    pw_diag_nmax: int = 50
    pw_diag_ndim: int = 4

    # === LCAO basis ===
    ks_solver: str = "genelpa"
    lcao_ecut: float = 0.0
    lcao_dk: float = 0.01
    lcao_dr: float = 0.01
    lcao_rmax: float = 30.0
    search_radius: float = -1.0
    search_pbc: bool = True
    gamma_only: int = 0

    # === Electronic structure ===
    nbands: int = 0
    nelec: float = 0.0
    dft_functional: str = "pbe"
    nspin: int = 1
    noncolin: bool = False
    lspinorb: bool = False
    smearing_method: str = "gauss"
    smearing_sigma: float = 0.01

    # === SCF ===
    scf_nmax: int = 100
    scf_thr: float = 1e-7
    scf_trust_min: float = 0.001
    scf_trust_max: float = 100.0
    mixing_type: str = "broyden"
    mixing_beta: float = 0.4
    mixing_ndim: int = 8
    mixing_gg0: float = 0.0
    chg_extrap: str = "atomic"

    # === Relaxation ===
    relax_nmax: int = 100
    relax_method: str = "cg"    # cg, bfgs, bfgs_trad, cg_bfgs, sd, fire
    relax_bfgs_ndim: int = 6
    relax_bfgs_wolfe: float = 0.01
    relax_bfgs_trust_radius_max: float = 0.8
    force_thr_ev: float = 0.01  # eV/Angstrom
    force_thr: float = 0.01     # eV/Angstrom (alias)
    stress_thr: float = 0.5     # kbar
    press1: float = 0.0
    press2: float = 0.0
    press3: float = 0.0
    fixed_axes: str = "None"
    cal_force: int = 1
    cal_stress: int = 1

    # === Molecular dynamics ===
    md_type: str = "nvt"        # nve, nvt, npt, f1, msst
    md_pmode: str = ""          # iso, aniso, tri (NPT pressure control mode, set only when md_type=npt)
    md_nstep: int = 1000
    md_dt: float = 1.0          # fs
    md_tfirst: float = 300.0    # K
    md_tlast: float = 300.0     # K
    md_damp: float = 0.5
    md_dumpfreq: int = 0        # MD_dump output frequency (steps, 0=unset)
    md_restartfreq: int = 0     # restart file output frequency (steps, 0=unset)
    out_level: str = ""         # output verbosity: ie (SCF) / i (relax) / m (MD simplified)
    dump_force: int = 0         # write forces to MD_dump
    dump_vel: int = 0           # write velocities to MD_dump

    # === Output control ===
    out_chg: bool = False
    out_pot: int = 0           # 0=off, 1=total local pot, 2=electrostatic pot
    out_dos: int = 0
    out_band: bool = False
    out_proj_band: bool = False
    out_bandgap: bool = False
    out_wfc_pw: bool = False
    out_wfc_lcao: bool = False
    out_wfc_r: bool = False
    out_elf: bool = False
    out_mat_hs: bool = False
    out_mat_hs2: bool = False
    out_dm: bool = False
    out_mul: bool = False
    out_pchg: bool = False

    # === DOS control ===
    dos_emin_ev: float = -15.0
    dos_emax_ev: float = 15.0
    dos_edelta_ev: float = 0.01
    dos_sigma: float = 0.07

    # === vdW ===
    vdw_method: str = "none"
    vdw_s6: float = 1.0
    vdw_s8: float = 1.0
    vdw_a1: float = 0.0
    vdw_a2: float = 0.0

    # === DFT+U ===
    dft_plus_u: bool = False
    orbital_corr: list[int] = field(default_factory=list)
    hubbard_u: list[float] = field(default_factory=list)

    # === Exact exchange (hybrid) ===
    exx_fock_alpha: float = 0.25
    exx_erfc_alpha: float = 0.0
    exx_erfc_omega: float = 0.0
    exx_separate_loop: bool = False
    exx_hybrid_step: int = 100

    # === TDDFT ===
    td_edm: int = 0
    td_propagator: int = 0
    td_stype: int = 0
    td_ttype: int = 0
    td_vext: bool = False

    # === Wannier ===
    berry_phase: bool = False
    gdir: int = 0
    towannier90: bool = False
    nnkpfile: str = ""

    # Catch-all for unrecognized parameters
    extras: dict[str, Any] = field(default_factory=dict, repr=False)

    def get_param(self, name: str, default=None):
        """Get a parameter value by name, checking both fields and extras."""
        if hasattr(self, name) and name not in ("extras", "_known_params"):
            return getattr(self, name)
        return self.extras.get(name, default)

    def set_param(self, name: str, value):
        """Set a parameter value, routing to field or extras."""
        if hasattr(self, name) and name not in ("extras", "_known_params"):
            setattr(self, name, value)
        else:
            self.extras[name] = value

    def as_dict(self) -> dict[str, Any]:
        """Return all parameters as a flat dict for writing to INPUT file."""
        result: dict[str, Any] = {}
        for field_name in self._known_params:
            val = getattr(self, field_name)
            # Skip fields with default/empty values to keep output clean
            if isinstance(val, bool):
                result[field_name] = val
            elif isinstance(val, list):
                if val:
                    result[field_name] = val
            elif val != 0 and val != 0.0 and val != "" and val is not None:
                result[field_name] = val
        result.update(self.extras)
        return result


# Auto-derived: every InputParams field name is a known parameter.
# Adding a new field to InputParams? It's automatically recognized — no manual sync needed.
InputParams._known_params = {
    f.name for f in fields(InputParams)
    if f.name not in ("extras", "_known_params")
}
