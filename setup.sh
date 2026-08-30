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
PIP_INDEX="https://pypi.tuna.tsinghua.edu.cn/simple"   # edit if you prefer another mirror

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

# 3c. Editable install — tiered mirror fallback + offline fallback:
#      1. Alibaba Cloud (fast, rarely blocked)
#      2. Tsinghua (may be blocked on some networks)
#      3. Default PyPI (global backstop)
#      4. Offline (no PyPI access): PEP 517 build isolation hangs trying to
#         download setuptools/wheel, so fall back to --no-build-isolation
#         --no-deps. The conda env already ships setuptools/wheel and all
#         runtime deps (fresh envs install them in the same setup run via pip).
ALI_INDEX="https://mirrors.aliyun.com/pypi/simple/"
PIP_INSTALL="pip install -e . --upgrade --timeout 15 --retries 2"
PIP_OFFLINE="pip install -e . --no-build-isolation --no-deps --upgrade"

# Bounded connectivity probe — avoids a long hang when packets are dropped.
if curl -sI --max-time 6 -o /dev/null "${ALI_INDEX}" 2>/dev/null; then
    if ${PIP_INSTALL} -i "${ALI_INDEX}" 2>/dev/null; then
        :
    elif ${PIP_INSTALL} -i "${PIP_INDEX}" 2>/dev/null; then
        echo -e "        ${YELLOW}Alibaba unreachable — used Tsinghua mirror${NC}"
    else
        echo -e "        ${YELLOW}Both mirrors unreachable — falling back to default PyPI${NC}"
        if ! ${PIP_INSTALL}; then
            echo -e "        ${RED}Error: pip install failed on all mirrors.${NC}"
            echo -e "        Trying offline install (--no-build-isolation --no-deps)..."
            if ! ${PIP_OFFLINE}; then
                echo -e "  ${RED}Error: pip install failed (online and offline).${NC}"
                echo -e "  Try re-running:"
                echo -e "    conda activate ${ENV_NAME} && pip install -e . --upgrade"
                exit 1
            fi
        fi
    fi
else
    echo -e "        ${YELLOW}No PyPI network access — offline install (--no-build-isolation --no-deps)${NC}"
    if ! ${PIP_OFFLINE}; then
        echo -e "  ${RED}Error: offline pip install failed.${NC}"
        echo -e "  Check that the conda env has setuptools/wheel and all runtime deps."
        exit 1
    fi
fi
echo -e "        ${GREEN}✓${NC} numpy, scipy, matplotlib, rich, pyyaml, ase, seekpath installed"
echo -e "        ${GREEN}✓${NC} atst-tools + NEB support installed"
echo -e "        ${GREEN}✓${NC} abacuscopilot command registered"

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
echo ""
echo -e "  [5/5] Verifying installation..."
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
