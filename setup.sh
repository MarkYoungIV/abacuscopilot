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

# 3c. Editable install — try the configured mirror first; fall back to
#     default PyPI if the mirror is unreachable (some machines can't
#     reach e.g. pypi.tuna.tsinghua.edu.cn).
ALI_INDEX="https://mirrors.aliyun.com/pypi/simple/"
PIP_INSTALL="pip install -e . --upgrade --timeout 15 --retries 2"

if ${PIP_INSTALL} -i "${PIP_INDEX}" 2>/dev/null; then
    :
elif ${PIP_INSTALL} -i "${ALI_INDEX}" 2>/dev/null; then
    echo -e "        ${YELLOW}Tsinghua unreachable — used Alibaba Cloud mirror${NC}"
else
    echo -e "        ${YELLOW}Both mirrors unreachable — falling back to default PyPI${NC}"
    if ! ${PIP_INSTALL}; then
        echo -e "  ${RED}Error: pip install failed on all mirrors.${NC}"
        echo -e "  Try re-running:"
        echo -e "    conda activate ${ENV_NAME} && pip install -e . --upgrade"
        exit 1
    fi
fi
echo -e "        ${GREEN}✓${NC} numpy, scipy, matplotlib, rich, pyyaml, ase, seekpath installed"
echo -e "        ${GREEN}✓${NC} atst-tools + NEB support installed"
echo -e "        ${GREEN}✓${NC} abacuscopilot command registered"

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
