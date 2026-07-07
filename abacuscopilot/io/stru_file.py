"""Reader and writer for ABACUS STRU files.

The STRU file defines the crystal structure: lattice parameters,
atomic positions, pseudopotential and orbital file names.

Format has five mandatory sections:
    ATOMIC_SPECIES
    NUMERICAL_ORBITAL  (omit for PW calculations)
    LATTICE_CONSTANT
    LATTICE_VECTORS
    ATOMIC_POSITIONS
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from abacuscopilot.core.exceptions import FileFormatError, FileNotFoundError_, MissingSectionError
from abacuscopilot.core.models import Atom, Structure


def read_stru(filepath: str | Path) -> Structure:
    """Read and parse an ABACUS STRU file.

    Args:
        filepath: Path to the STRU file.

    Returns:
        Structure object with all parsed data.

    Raises:
        FileNotFoundError_: If file doesn't exist.
        FileFormatError: If the file format is invalid.
        MissingSectionError: If a required section is missing.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError_(str(filepath), "STRU file not found")

    with open(filepath) as f:
        lines = f.readlines()

    # Clean lines: strip comments and whitespace
    clean_lines = []
    for line in lines:
        # Remove inline comments (everything after # or !)
        for comment_char in ("#", "!"):
            idx = line.find(comment_char)
            if idx >= 0:
                line = line[:idx]
        line = line.strip()
        if line:
            clean_lines.append(line)

    return _parse_stru_lines(clean_lines, str(filepath))


def _parse_stru_lines(lines: list[str], filepath: str) -> Structure:
    """Parse cleaned STRU lines into a Structure object."""
    structure = Structure()

    # Find section boundaries
    idx = 0
    section_order = []

    while idx < len(lines):
        line = lines[idx]
        # Strip leading line numbers from copy-paste artifacts (e.g. "8 LATTICE_CONSTANT")
        clean = re.sub(r"^\d+\s+", "", line)
        upper = clean.upper()

        if upper.startswith("ATOMIC_SPECIES"):
            section_order.append(("ATOMIC_SPECIES", idx + 1))
        elif upper.startswith("NUMERICAL_ORBITAL"):
            section_order.append(("NUMERICAL_ORBITAL", idx + 1))
        elif upper.startswith("LATTICE_CONSTANT"):
            section_order.append(("LATTICE_CONSTANT", idx + 1))
        elif upper.startswith("LATTICE_VECTORS"):
            section_order.append(("LATTICE_VECTORS", idx + 1))
        elif upper.startswith("ATOMIC_POSITIONS"):
            section_order.append(("ATOMIC_POSITIONS", idx + 1))

        idx += 1

    # Parse each section
    section_data = {}
    for i, (sec_name, start) in enumerate(section_order):
        end = section_order[i + 1][1] - 1 if i + 1 < len(section_order) else len(lines)
        section_data[sec_name] = lines[start:end]

    # --- ATOMIC_SPECIES ---
    if "ATOMIC_SPECIES" not in section_data:
        raise MissingSectionError(filepath, "ATOMIC_SPECIES")

    for line in section_data["ATOMIC_SPECIES"]:
        parts = line.split()
        if len(parts) >= 2:
            label = parts[0]
            if label not in structure.species_order:
                structure.species_order.append(label)
            if len(parts) >= 3:
                structure.pseudo_files[label] = parts[2]

    # --- NUMERICAL_ORBITAL (optional) ---
    if "NUMERICAL_ORBITAL" in section_data:
        lines_orb = section_data["NUMERICAL_ORBITAL"]
        for i, line in enumerate(lines_orb):
            if i < len(structure.species_order):
                structure.orbital_files[structure.species_order[i]] = line
            else:
                structure.orbital_files[f"species_{i}"] = line

    # --- LATTICE_CONSTANT ---
    if "LATTICE_CONSTANT" not in section_data:
        raise MissingSectionError(filepath, "LATTICE_CONSTANT")

    const_line = section_data["LATTICE_CONSTANT"][0]
    structure.lattice.constant = float(const_line.split()[0])

    # --- LATTICE_VECTORS ---
    if "LATTICE_VECTORS" not in section_data:
        raise MissingSectionError(filepath, "LATTICE_VECTORS")

    vectors = np.zeros((3, 3))
    for i, line in enumerate(section_data["LATTICE_VECTORS"][:3]):
        parts = line.split()
        if len(parts) >= 3:
            vectors[i] = [float(x) for x in parts[:3]]
    structure.lattice.vectors = vectors

    # --- ATOMIC_POSITIONS ---
    if "ATOMIC_POSITIONS" not in section_data:
        raise MissingSectionError(filepath, "ATOMIC_POSITIONS")

    pos_lines = section_data["ATOMIC_POSITIONS"]
    if not pos_lines:
        raise FileFormatError(filepath, "ATOMIC_POSITIONS section is empty")

    # First line: coordinate type
    structure.coordinate_type = pos_lines[0]

    pi = 1
    while pi < len(pos_lines):
        # Element label line
        element = pos_lines[pi]
        pi += 1

        if pi >= len(pos_lines):
            break

        # Default magnetism
        mag_line = pos_lines[pi]
        try:
            default_mag = float(mag_line)
            structure.magnetism[element] = default_mag
        except ValueError:
            default_mag = 0.0
            pi -= 1  # This wasn't a magnetism line, rewind
        pi += 1

        if pi >= len(pos_lines):
            break

        # Number of atoms
        try:
            num_atoms = int(pos_lines[pi])
        except ValueError as e:
            raise FileFormatError(
                filepath,
                f"Expected number of atoms for element {element}, got '{pos_lines[pi]}'"
            ) from e
        pi += 1

        # Atom positions
        for _ in range(num_atoms):
            if pi >= len(pos_lines):
                raise FileFormatError(
                    filepath,
                    f"Expected {num_atoms} positions for {element}, but reached end of section"
                )

            parts = pos_lines[pi].split()
            if len(parts) < 3:
                raise FileFormatError(
                    filepath,
                    f"Atom position must have at least 3 coordinates, got: {pos_lines[pi]}"
                )

            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])

            # Parse optional keywords after coordinates
            fix = (True, True, True)
            magmom = 0.0
            velocity = None
            angle1 = None
            angle2 = None

            ki = 3
            while ki < len(parts):
                kw = parts[ki].lower()
                if kw in ("m",):
                    # Movement flags: m v1 v2 v3 (0=fixed, 1=free)
                    flags = []
                    for _ in range(3):
                        ki += 1
                        if ki < len(parts):
                            flags.append(parts[ki] not in ("0",))
                    if len(flags) == 3:
                        fix = tuple(flags)
                elif kw in ("mag", "magmom"):
                    ki += 1
                    if ki < len(parts):
                        magmom = float(parts[ki])
                elif kw == "v":
                    vel = []
                    for _ in range(3):
                        ki += 1
                        if ki < len(parts):
                            vel.append(float(parts[ki]))
                    if len(vel) == 3:
                        velocity = np.array(vel)
                elif kw == "angle1":
                    ki += 1
                    if ki < len(parts):
                        angle1 = float(parts[ki])
                elif kw == "angle2":
                    ki += 1
                    if ki < len(parts):
                        angle2 = float(parts[ki])
                ki += 1

            structure.atoms.append(
                Atom(
                    species=element,
                    position=np.array([x, y, z]),
                    fix=fix,
                    magmom=magmom,
                    velocity=velocity,
                    angle1=angle1,
                    angle2=angle2,
                )
            )
            pi += 1

    return structure


# Atomic masses (g/mol) — standard IUPAC values
_ATOMIC_MASSES: dict[str, float] = {
    "H": 1.008, "He": 4.002602, "Li": 6.941, "Be": 9.012183, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.1797,
    "Na": 22.989769, "Mg": 24.305, "Al": 26.9815385, "Si": 28.085,
    "P": 30.973762, "S": 32.065, "Cl": 35.453, "Ar": 39.948,
    "K": 39.0983, "Ca": 40.078, "Sc": 44.955908, "Ti": 47.867,
    "V": 50.9415, "Cr": 51.9961, "Mn": 54.938043, "Fe": 55.845,
    "Co": 58.933194, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.63, "As": 74.921595, "Se": 78.971,
    "Br": 79.904, "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62,
    "Y": 88.905838, "Zr": 91.224, "Nb": 92.90637, "Mo": 95.95,
    "Tc": 98.0, "Ru": 101.07, "Rh": 102.90549, "Pd": 106.42,
    "Ag": 107.8682, "Cd": 112.414, "In": 114.818, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.6, "I": 126.90447, "Xe": 131.293,
    "Cs": 132.905452, "Ba": 137.327, "La": 138.90547, "Ce": 140.116,
    "Pr": 140.90766, "Nd": 144.242, "Pm": 145.0, "Sm": 150.36,
    "Eu": 151.964, "Gd": 157.25, "Tb": 158.925354, "Dy": 162.5,
    "Ho": 164.930328, "Er": 167.259, "Tm": 168.934218, "Yb": 173.054,
    "Lu": 174.9668, "Hf": 178.49, "Ta": 180.94788, "W": 183.84,
    "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
    "Au": 196.96657, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2,
    "Bi": 208.9804, "Po": 209.0, "At": 210.0, "Rn": 222.0,
    "Fr": 223.0, "Ra": 226.0, "Ac": 227.0, "Th": 232.0377,
    "Pa": 231.03588, "U": 238.02891, "Np": 237.0, "Pu": 244.0,
}


def write_stru(
    structure: Structure,
    filepath: str | Path = "STRU",
    is_lcao: bool = True,
    is_dp: bool = False,
) -> None:
    """Write a Structure object to an ABACUS STRU file.

    Args:
        is_dp: If True, omit pseudopotential and orbital filenames
               (Deep Potential mode, no DFT needed).

    Output format matches the abacustest convention:
    - Atomic masses from IUPAC standard values
    - LATTICE_CONSTANT: no indent, 6 decimal places
    - LATTICE_VECTORS: 4-space indent, 11 decimal places
    - ATOMIC_POSITIONS: 4-space indent for coordinates, movement flags (1 1 1)

    Args:
        structure: Structure object to write.
        filepath: Output file path.
        is_lcao: Whether to include NUMERICAL_ORBITAL section (LCAO only).
    """
    filepath = Path(filepath)
    lines = []

    # === ATOMIC_SPECIES ===
    lines.append("ATOMIC_SPECIES")
    for species in structure.species_order:
        mass = _ATOMIC_MASSES.get(species, 0.0)
        if is_dp:
            lines.append(f"{species} {mass:.6f}")
        else:
            pseudo = structure.pseudo_files.get(species, f"{species}.upf")
            lines.append(f"{species} {mass:.6f} {pseudo}")
    lines.append("")

    # === NUMERICAL_ORBITAL (LCAO only; skip for DP and PW) ===
    if is_lcao and not is_dp:
        lines.append("NUMERICAL_ORBITAL")
        for species in structure.species_order:
            orb = structure.orbital_files.get(species, f"{species}.orb")
            lines.append(f"{orb}")
        lines.append("")

    # === LATTICE_CONSTANT ===
    lines.append("LATTICE_CONSTANT")
    lines.append(f"{structure.lattice.constant:.6f}")
    lines.append("")

    # === LATTICE_VECTORS ===
    lines.append("LATTICE_VECTORS")
    for row in structure.lattice.vectors:
        lines.append(f"    {row[0]:15.10f}  {row[1]:15.10f}  {row[2]:15.10f}")
    lines.append("")

    # === ATOMIC_POSITIONS ===
    lines.append("ATOMIC_POSITIONS")
    lines.append(f"{structure.coordinate_type}")
    lines.append("")

    for species in structure.species_order:
        atoms_of_species = [a for a in structure.atoms if a.species == species]
        if not atoms_of_species:
            continue

        # Element label
        lines.append(species)

        # Default magnetism
        default_mag = structure.magnetism.get(species, 0.0)
        lines.append(f"{default_mag:.6f}")

        # Number of atoms
        lines.append(f"{len(atoms_of_species)}")

        # Positions — fixed-width for alignment regardless of magnitude
        for atom in atoms_of_species:
            x, y, z = atom.position[0], atom.position[1], atom.position[2]
            pos_str = f"    {x:15.10f}  {y:15.10f}  {z:15.10f}"

            # Movement flags with ABACUS standard 'm' prefix
            flags = " ".join("1" if f else "0" for f in atom.fix)
            pos_str += f" m {flags}"

            # Magnetic moment
            if atom.magmom != 0.0:
                pos_str += f" mag {atom.magmom:.6f}"

            # Velocity
            if atom.velocity is not None:
                v = atom.velocity
                pos_str += f" v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}"

            # Non-collinear angles
            if atom.angle1 is not None:
                pos_str += f" angle1 {atom.angle1:.6f}"
            if atom.angle2 is not None:
                pos_str += f" angle2 {atom.angle2:.6f}"

            lines.append(pos_str)
        lines.append("")

    with open(filepath, "w") as f:
        f.write("\n".join(lines))
