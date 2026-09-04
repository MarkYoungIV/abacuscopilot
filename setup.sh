#!/usr/bin/env bash
# =============================================================================
# AbacusCopilot Setup Script
# =============================================================================
# One-command install / upgrade for macOS and Linux.
#
# Idempotent: safe to run whether or not abacuscopilot was installed before.
#   - Creates the conda env if missing, reuses it if present.
#   - On upgrade, removes stale pip registration + bytecode caches so files
#     deleted in a new version do not linger as ghosts.
#
# Usage:
#   chmod +x setup.sh
#   ./setup.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="abacuscopilot"
PY_VERSION="3.11"
# Public PyPI mirrors are configured in the MIRRORS list under section 3c.

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

echo ""
echo -e "  ${CYAN}${BOLD}AbacusCopilot Setup${NC}"
echo -e "  ${CYAN}A pre- & post-processing toolkit for ABACUS DFT${NC}"
echo ""

# --- 0. Preflight: conda must be available ---
if ! command -v conda >/dev/null 2>&1; then
    echo -e "  ${RED}Error: 'conda' not found in PATH.${NC}"
    echo -e "  Please install Anaconda or Miniconda first, then re-run this script."
    exit 1
fi

# --- 1. Detect OS ---
OS_NAME="Linux"
if [[ "$(uname)" == "Darwin" ]]; then
    OS_NAME="macOS"
fi
echo -e "  [1/5] Detected OS: ${GREEN}${OS_NAME}${NC}"

# --- 2. Setup conda environment (reuse if present, create if missing) ---
# Exact whole-word match on the env name column, so 'abacuscopilot_dev' etc.
# never gets mistaken for 'abacuscopilot'.
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo -e "  [2/5] Conda env '${ENV_NAME}' already exists — reusing it"
    IS_UPGRADE=1
else
    echo -e "  [2/5] Creating conda env '${ENV_NAME}' with Python ${PY_VERSION}..."
    conda create -n "${ENV_NAME}" python="${PY_VERSION}" -y
    IS_UPGRADE=0
fi

# Activate the env.
# Temporarily disable `set -u` (nounset): conda's own activate/deactivate.d
# hook scripts (e.g. from gxx_linux-64) reference unset vars like
# CONDA_BACKUP_CXX, which would abort the script under `set -u`.
set +u
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME}"
set -u
PY_VER=$(python --version 2>&1)
echo -e "        ${GREEN}${PY_VER}${NC}"

# --- 3. Install / upgrade abacuscopilot ---
echo ""
echo -e "  [3/5] Installing abacuscopilot + dependencies..."
cd "$SCRIPT_DIR"

# 3a. Clear stale bytecode caches in the source tree (avoid ghost .pyc from an
#     old version whose .py files were removed).
find "$SCRIPT_DIR" -name '__pycache__' -type d -prune -exec rm -rf {} + 2>/dev/null || true
find "$SCRIPT_DIR" -name '*.pyc' -delete 2>/dev/null || true

# 3b. On upgrade, uninstall the previously registered package first so pip's
#     record of old files is dropped cleanly before re-installing.
if [[ "${IS_UPGRADE}" == "1" ]] && pip show abacuscopilot >/dev/null 2>&1; then
    echo -e "        ${YELLOW}Existing install detected — removing it first${NC}"
    pip uninstall -y abacuscopilot >/dev/null 2>&1 || true
fi

# 3c. Editable install — try, in order:
#   1. User's own pip index (pip.conf / PIP_INDEX_URL) — internal-network
#      machines (e.g. HPC compute nodes) usually have a mirror configured that
#      no public list knows about, so honor it first.
#   2. Public PyPI mirrors — probe + install with the first one that both
#      responds and actually installs cleanly.
#   3. Offline fallback (no PyPI access): PEP 517 build isolation hangs trying
#      to download setuptools/wheel, so fall back to --no-build-isolation
#      --no-deps. The conda env ships setuptools/wheel; runtime deps are then
#      verified below and, if missing, fail loudly with guidance.
PIP_INSTALL="pip install -e . --upgrade --timeout 15 --retries 2"
PIP_OFFLINE="pip install -e . --no-build-isolation --no-deps --upgrade"

# Public PyPI mirrors, ordered by reachability for China mainland.
# Edit this list if you prefer different mirrors.
MIRRORS=(
    "https://mirrors.aliyun.com/pypi/simple/"
    "https://pypi.tuna.tsinghua.edu.cn/simple"
    "https://mirrors.ustc.edu.cn/pypi/simple/"
    "https://mirrors.huaweicloud.com/repository/pypi/simple/"
    "https://mirrors.cloud.tencent.com/pypi/simple/"
    "https://pypi.org/simple"
)

DEPS_INSTALLED=0

# 1. Honor the user's existing pip index first (covers institutional mirrors).
USER_INDEX="${PIP_INDEX_URL:-}"
[ -z "$USER_INDEX" ] && USER_INDEX="$(pip config get global.index-url 2>/dev/null | tr -d '[:space:]' || true)"
[ "$USER_INDEX" = "None" ] && USER_INDEX=""
if [ -n "$USER_INDEX" ] && curl -sI --max-time 6 -o /dev/null "$USER_INDEX" 2>/dev/null; then
    echo -e "        Using configured pip index: ${CYAN}${USER_INDEX}${NC}"
    if ${PIP_INSTALL} -i "$USER_INDEX" 2>/dev/null; then
        echo -e "        ${GREEN}✓${NC} Installed from your configured index"
        DEPS_INSTALLED=1
    else
        echo -e "        ${YELLOW}Install failed with your configured index — trying public mirrors${NC}"
    fi
fi

# 2. Otherwise (or on failure above) probe + install with each public mirror.
if [ "${DEPS_INSTALLED}" != "1" ]; then
    for M in "${MIRRORS[@]}"; do
        if curl -sI --max-time 6 -o /dev/null "$M" 2>/dev/null && ${PIP_INSTALL} -i "$M" 2>/dev/null; then
            echo -e "        ${GREEN}✓${NC} Installed from mirror: ${CYAN}${M}${NC}"
            DEPS_INSTALLED=1
            break
        fi
        echo -e "        ${YELLOW}✗ ${M} unreachable or install failed${NC}"
    done
fi

# 3. All online paths failed — offline install, then verify deps honestly.
if [ "${DEPS_INSTALLED}" != "1" ]; then
    echo -e "  ${YELLOW}No reachable PyPI mirror — trying offline install (--no-build-isolation --no-deps)${NC}"
    if ${PIP_OFFLINE}; then
        # --no-deps skips every runtime dependency; check and fail loudly rather
        # than finishing with a broken install (this was silent before).
        if python -c "import numpy, scipy, matplotlib, rich, yaml, ase, seekpath, phonopy" 2>/dev/null; then
            DEPS_INSTALLED=1
        else
            echo -e "  ${RED}ERROR: installed without dependencies (offline mode).${NC}"
            echo -e "  This server cannot reach any PyPI mirror. Install deps by:"
            echo -e "    1) On a networked machine: pip download <deps> -d wheels/, copy over, then"
            echo -e "       pip install wheels/*.whl && ${PIP_OFFLINE}"
            echo -e "    2) Or via conda-forge: conda install -c conda-forge numpy scipy matplotlib ase phonopy"
            exit 1
        fi
    else
        echo -e "  ${RED}Error: pip install failed (online and offline).${NC}"
        echo -e "  Try re-running:"
        echo -e "    conda activate ${ENV_NAME} && pip install -e . --upgrade"
        exit 1
    fi
fi
echo -e "        ${GREEN}✓${NC} abacuscopilot + dependencies installed"

# 3d. Build the bundled Bader program (needed by Bader Charge, task 1304;
#     currently hidden pending validation).
#     Prefer compiling with gfortran (matches the target arch); if gfortran is
#     missing or the compile fails, download the prebuilt binary for the
#     platform (Linux x86-64 / macOS) from the Henkelman site.
if [ -d "$SCRIPT_DIR/scripts/bader/bader" ]; then
    _BADER_X="$SCRIPT_DIR/scripts/bader/bader.x"
    echo -e "        Building Bader charge program..."
    if command -v gfortran >/dev/null 2>&1; then
        if ( cd "$SCRIPT_DIR/scripts/bader/bader" && make -f makefile.osx_gfortran 2>/dev/null ) || \
           ( cd "$SCRIPT_DIR/scripts/bader/bader" && make -f makefile.osx_gfortran LINK="" 2>/dev/null ); then
            cp -f "$SCRIPT_DIR/scripts/bader/bader/bader" "$_BADER_X"
            chmod +x "$_BADER_X"
        fi
    fi
    # Wrong-arch or missing binary (e.g. a shipped macOS arm64 binary on a
    # Linux box, or no gfortran) → download the prebuilt for this platform.
    if [ ! -x "$_BADER_X" ] || ! $_BADER_X 2>&1 | grep -q "BADER"; then
        case "$(uname -s)" in
            Linux)  _BADER_PKG="bader_lnx_64.tar.gz" ;;
            Darwin) _BADER_PKG="bader_osx_gfortran.tar.gz" ;;
            *)      _BADER_PKG="" ;;
        esac
        if [ -n "$_BADER_PKG" ]; then
            echo -e "        Downloading prebuilt Bader ($_BADER_PKG)..."
            _BADER_URL="https://theory.cm.utexas.edu/henkelman/code/bader/download/$_BADER_PKG"
            _tmpd="$(mktemp -d)"
            if curl -fsSL --max-time 120 -o "$_tmpd/$_BADER_PKG" "$_BADER_URL" \
               && ( cd "$_tmpd" && tar xzf "$_BADER_PKG" ) \
               && cp -f "$_tmpd/bader" "$_BADER_X" && chmod +x "$_BADER_X"; then
                echo -e "        ${GREEN}✓${NC} Bader program ready (scripts/bader/bader.x)"
            else
                echo -e "        ${YELLOW}! Could not build/download Bader — Bader Charge (1304) will be unavailable.${NC}"
            fi
            rm -rf "$_tmpd"
        fi
    else
        echo -e "        ${GREEN}✓${NC} Bader program ready (scripts/bader/bader.x)"
    fi
fi

# --- 4. Generate default config ---
echo ""
echo -e "  [4/5] Configuring..."
if [ ! -f "$HOME/.abacuscopilot/config.yaml" ]; then
    python -c "from abacuscopilot.config import load_config; load_config()"
    echo -e "        ${GREEN}✓${NC} Created ~/.abacuscopilot/config.yaml"
    echo -e "        ${GREEN}✓${NC} Pseudo/orbital paths auto-detected from package"
else
    echo -e "        ${GREEN}✓${NC} Config already exists (not overwritten)"
fi

# --- 4b. Library availability ---
# Pseudopotential/orbital libraries (PP-Orb/) are NOT bundled in the git repo
# (large + third-party redistribution licensing). A fresh clone has none, while
# the full release tarball ships them. Hint the user instead of failing silently.
if [ ! -d "$SCRIPT_DIR/PP-Orb" ] || [ -z "$(ls -A "$SCRIPT_DIR/PP-Orb" 2>/dev/null)" ]; then
    echo ""
    echo -e "  ${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "  ${YELLOW}  Pseudopotential/orbital libraries (PP-Orb/) not found.${NC}"
    echo -e "  ${YELLOW}  They are NOT bundled in the git repo. Options:${NC}"
    echo -e "  ${YELLOW}    1. Download the full release tarball (abacuscopilot_v*.tar.gz)${NC}"
    echo -e "  ${YELLOW}       from the GitHub Releases page — it contains PP-Orb/ — and${NC}"
    echo -e "  ${YELLOW}       re-run ./setup.sh from that extracted directory; or${NC}"
    echo -e "  ${YELLOW}    2. Place your SG15 / lanthanide UPF + orbital files under${NC}"
    echo -e "  ${YELLOW}       PP-Orb/ yourself.${NC}"
    echo -e "  ${YELLOW}  AbacusCopilot still installs and runs without them, but warns${NC}"
    echo -e "  ${YELLOW}  when generating INPUT/STRU for elements it cannot resolve.${NC}"
    echo -e "  ${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
fi

# --- 5. Verify installed version ---
# Fail loudly if the import breaks, instead of printing "version: unknown"
# and "Setup complete!" over a broken install.
echo ""
echo -e "  [5/5] Verifying installation..."
if ! python -c "import abacuscopilot" >/dev/null 2>&1; then
    echo -e "  ${RED}ERROR: abacuscopilot import failed — install is incomplete.${NC}"
    echo -e "  Run:  conda activate ${ENV_NAME} && pip install -e . --no-build-isolation"
    exit 1
fi
INSTALLED_VER=$(python -c "import abacuscopilot; print(abacuscopilot.__version__)" 2>/dev/null || echo "unknown")
echo -e "        ${GREEN}✓${NC} AbacusCopilot version: ${BOLD}${INSTALLED_VER}${NC}"

echo ""
echo -e "  ${GREEN}${BOLD}Setup complete!${NC}"
echo ""
echo -e "  To start:  ${BOLD}conda activate ${ENV_NAME} && abacuscopilot${NC}"
echo ""
echo -e "  ${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "  ${YELLOW}  MLP-SSCHA (tasks 3205-3208) requires extra packages:${NC}"
echo -e "  ${YELLOW}    pypolymlp  symfc  dpdata${NC}"
echo -e "  ${YELLOW}  Recommended install (with Tsinghua mirror):${NC}"
echo -e "  ${YELLOW}    conda install pypolymlp symfc -c conda-forge${NC}"
echo -e "  ${YELLOW}    pip install dpdata -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple${NC}"
echo -e "  ${YELLOW}  If conda is slow, add mirror first:${NC}"
echo -e "  ${YELLOW}    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge${NC}"
echo -e "  ${YELLOW}  These are NOT required for standard phonon (112/113) or other tasks.${NC}"
echo -e "  ${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
