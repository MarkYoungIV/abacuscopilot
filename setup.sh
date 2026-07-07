#!/usr/bin/env bash
# =============================================================================
# AbacusCopilot Setup Script
# =============================================================================
# One-command installation for macOS and Linux.
#
# Usage:
#   chmod +x setup.sh
#   ./setup.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

echo ""
echo -e "  ${CYAN}${BOLD}AbacusCopilot Setup${NC}"
echo -e "  ${CYAN}A pre- & post-processing toolkit for ABACUS DFT${NC}"
echo ""

# --- 1. Detect OS ---
OS_NAME="Linux"
if [[ "$(uname)" == "Darwin" ]]; then
    OS_NAME="macOS"
fi
echo -e "  [1/4] Detected OS: ${GREEN}${OS_NAME}${NC}"

# --- 2. Setup conda environment ---
ENV_NAME="abacuscopilot"
if conda env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
    echo -e "  [2/4] Conda env '${ENV_NAME}' already exists — using it"
else
    echo -e "  [2/4] Creating conda env '${ENV_NAME}' with Python 3.11..."
    conda create -n ${ENV_NAME} python=3.11 -y
fi

eval "$(conda shell.bash hook)"
conda activate ${ENV_NAME}
PY_VER=$(python --version 2>&1)
echo -e "        ${GREEN}${PY_VER}${NC}"

# --- 3. Install abacuscopilot ---
echo ""
echo -e "  [3/4] Installing abacuscopilot + dependencies..."
cd "$SCRIPT_DIR"
pip install -e . -i https://pypi.tuna.tsinghua.edu.cn/simple --quiet 2>&1 | tail -1
echo -e "        ${GREEN}✓${NC} numpy, scipy, matplotlib, rich, pyyaml installed"
echo -e "        ${GREEN}✓${NC} abacuscopilot command registered"

# --- 4. Generate default config ---
echo ""
echo -e "  [4/4] Configuring..."
if [ ! -f "$HOME/.abacuscopilot/config.yaml" ]; then
    python -c "from abacuscopilot.config import load_config; load_config()"
    echo -e "        ${GREEN}✓${NC} Created ~/.abacuscopilot/config.yaml"
    echo -e "        ${GREEN}✓${NC} Pseudo/orbital paths auto-detected from package"
else
    echo -e "        ${GREEN}✓${NC} Config already exists (not overwritten)"
fi

echo ""
echo -e "  ${GREEN}${BOLD}Setup complete!${NC}"
echo ""
echo -e "  To start:  ${BOLD}conda activate ${ENV_NAME} && abacuscopilot${NC}"
echo ""
