# Changelog

All notable changes to AbacusCopilot are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.38.4] - 2026-09-23

### Fixed

- **Lattice redefinition (task 402) applied the transformation matrix on the
  wrong side.** Lattice vectors are the *rows* of L, so the new cell is
  `M @ L` — row i of M gives the combination of old vectors forming the i-th
  new vector (a' = m₁₁·a + m₁₂·b + m₁₃·c). The code used `L @ M`, which for an
  orthogonal cell merely rescales each row of M by the corresponding lattice
  constant. The resulting cell vectors were therefore *not* integer
  combinations of the old ones (coefficients like 0.858 / −1.165), so
  `Redefined.STRU` was not a supercell of the input crystal at all. The
  fractional-coordinate transform was transposed to match
  (`f' = f @ M⁻¹`). Volume ratio and atom count were correct either way, which
  is why the bug stayed silent — VASPKIT, given the same matrix, produced a
  different (correct) cell for the same input. The convention now matches
  VASPKIT's TRANSMAT, so integer matrices can be copied across verbatim.

- **Supercell construction (task 401) had the same defect** and is fixed to
  `scale @ L`. For an orthogonal cell the two orders coincide, so only
  non-orthogonal (monoclinic / trigonal / triclinic) input was affected: the
  cell came out with the Cartesian components scaled rather than the lattice
  vectors.

- The replica-offset search in task 402 now bounds its search box by the *row*
  sums of |M| (allowing negative components) instead of the column sums over a
  non-negative range. Previously a matrix with negative entries could yield
  fewer than |det M| offsets, silently disabling atom replication and writing
  a cell with too few atoms.

### Added

- Regression tests for both tasks (`tests/test_structure_editing.py`): the new
  vectors must be integer combinations of the old, fractional coordinates must
  match VASPKIT's `SUPERCELL.vasp` for the same input, and every replica must
  land on an original lattice site — checked on a non-orthogonal cell too,
  where an orthogonal-only shortcut cannot pass.

## [0.1.38c] - 2026-09-22

### Fixed

- **Resource resolution now defaults to SG15 in every flow.**  The active
  `libraries.family` owns the library lists: with the default `sg15`,
  `pseudo_library` / `orbital_library` are re-derived from the bundled
  SG15 + lanthanide roots on every load, so stale entries from another series
  (e.g. a Dojo dir left in the list) can no longer supply files — a Si CIF
  used to resolve to Dojo's `Si.upf` + SZ `Si_gga_7au_100Ry_1s1p.orb` instead
  of SG15's `Si_ONCV_PBE-1.0.upf` + `2s2p1d`.  `custom` keeps hand-configured
  lists; an explicitly picked family (APNS, Dojo) is unaffected.

- **The default-basis tie-break no longer picks the smallest orbital.**
  `_candidate_rank` keyed on one hardcoded DZP name (`4s2p2d1f`), which no
  main-group file matches; every Si candidate tied and the alphabetical
  fallback chose SZ `1s1p` over DZP `2s2p1d`.  Ranking now prefers the
  canonical `4s2p2d1f`, then any parseable DZP-like basis, then
  unparseable/minimal names, then rcut closest to 7 au, then name.

### Added

- Regression tests for family-owned library lists (`tests/test_config.py`)
  and the DZP-before-SZ ranking (`tests/test_library_resolution.py`).

## [0.1.38b] - 2026-09-21

### Fixed

- **Elastic constants (task 1201): general Born stability criterion.**
  The stability verdict no longer applies the cubic inequalities
  (C₁₁−C₁₂>0, C₁₁+2C₁₂>0, C₄₄>0) to every crystal.  The relaxed-ion 6×6
  elastic matrix is now tested for positive definiteness (all eigenvalues
  > 0) — the necessary-and-sufficient criterion, valid for any crystal
  system (Mouhat & Coudert, Phys. Rev. B 90, 224104 (2014)).  All six
  eigenvalues and the softest strain mode are printed, and
  `elastic_properties.dat` gains a minimum-eigenvalue line.

- **EOS Setup (task 110) now generates `calculation relax` INPUT files.**
  Scaling the lattice while keeping the original fractional coordinates in a
  plain SCF run gives a too-stiff E(V) curve and an overestimated B₀.  Each
  `scale_*` directory now relaxes the ions at the fixed scaled cell, and EOS
  Fitting (task 1202) reads `OUT.ABACUS/running_relax.log`, falling back to
  `running_scf.log` for existing SCF-based directories.

## [0.1.38] - 2026-09-20

### Fixed

- **MD trajectory analysis — triclinic (NPT) cell handling (task 3107).**
  The LAMMPS dump → ABACUS MD_dump conversion now parses the `BOX BOUNDS`
  tilt columns (`xy`, `xz`, `yz`), converts the LAMMPS *bounding-box* values
  to the true box edges, and writes the full 3×3 lattice. Previously the tilt
  columns were dropped and the bounding-box lengths were written as an
  orthogonal (diagonal) cell, which overestimated the cell volume by
  10–20% for sheared boxes and discarded the cell shape.

- **Minimum-image convention in RDF and pairwise distances (task 3104).**
  `_minimum_image` (and the pairwise-distance helper) used the transposed
  form `inv(cell.T)` / `@ cell.T`, whereas the code base consistently uses the
  row-vector convention `cart = frac @ cell`. Wrapping therefore happened in
  the wrong space for any non-orthogonal cell. The implementation now uses
  `inv(cell)` / `@ cell` and, for triclinic cells, searches the 26 neighbouring
  images for the true shortest image (matching ASE's
  `get_all_distances(mic=True)` and a brute-force nearest-image search). The
  PyTorch/MPS fast path is now restricted to orthogonal cells, where
  box-length wrapping is exact.

- **Same-convention fixes in the remaining MD/PBC helpers.** The corrected
  row-vector minimum image is also applied to POSCAR export (`write_poscar`),
  probability-density gridding, MSD trajectory unwrapping
  (`_unwrap_trajectory`), and the van Hove non-Gaussian parameter.

### Added

- Regression tests for triclinic-cell handling (`tests/test_md_triclinic.py`):
  minimum image vs. brute force (mild and highly sheared cells), LAMMPS
  `BOX BOUNDS` → cell conversion (triclinic, orthorhombic, positive tilt),
  pairwise distances, and an end-to-end task 3107 dump conversion.

## [0.1.37] - 2026-09-12

- Previous release.
