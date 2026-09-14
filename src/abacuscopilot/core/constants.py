"""Physical constants and conversion factors used across abacuscopilot.

All values are in atomic units (Hartree) or derived SI unless noted otherwise.
"""

# === Energy conversions ===
RY_TO_EV = 13.605693122994    # Rydberg to eV
EV_TO_RY = 1.0 / RY_TO_EV     # eV to Rydberg
HA_TO_EV = 27.211386245988    # Hartree to eV
EV_TO_HA = 1.0 / HA_TO_EV     # eV to Hartree
RY_TO_HA = 0.5                # Rydberg to Hartree
HA_TO_RY = 2.0                # Hartree to Rydberg

# === Length conversions ===
BOHR_TO_ANGSTROM = 0.529177210903   # Bohr to Angstrom
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM  # Angstrom to Bohr

# === Force conversions ===
RY_PER_BOHR_TO_EV_PER_ANGSTROM = RY_TO_EV / BOHR_TO_ANGSTROM
EV_PER_ANGSTROM_TO_RY_PER_BOHR = 1.0 / RY_PER_BOHR_TO_EV_PER_ANGSTROM

# === Pressure conversions ===
RY_PER_BOHR3_TO_KBAR = 14710.507848350711  # Ry/Bohr^3 to kbar
KBAR_TO_RY_PER_BOHR3 = 1.0 / RY_PER_BOHR3_TO_KBAR
KBAR_TO_GPA = 0.1                          # kbar to GPa
GPA_TO_KBAR = 10.0                         # GPa to kbar

# === Mass conversions ===
AMU_TO_KG = 1.66053906660e-27             # atomic mass unit to kg
ELECTRON_MASS_KG = 9.1093837015e-31       # electron rest mass in kg

# === Time conversions ===
AU_TIME_TO_FS = 0.02418884254              # atomic time unit to femtoseconds
FS_TO_AU_TIME = 1.0 / AU_TIME_TO_FS

# === Temperature ===
KB_EV = 8.617333262145e-5                 # Boltzmann constant in eV/K

# === Common ABACUS defaults ===
DEFAULT_ECUTWFC = 100.0                    # Ry
DEFAULT_KSPACING = 0.14                    # 1/bohr — ABACUS kspacing unit (manual: suggest < 0.25)
DEFAULT_SCF_THR = 1e-7                     # Ry
DEFAULT_FORCE_THR = 0.001                  # eV/Å
DEFAULT_STRESS_THR = 0.5                   # kbar
DEFAULT_SCF_NMAX = 100
DEFAULT_RELAX_NMAX = 100
DEFAULT_MD_DT = 1.0                        # fs
DEFAULT_MD_NSTEP = 1000
