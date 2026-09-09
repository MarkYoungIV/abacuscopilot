#!/usr/bin/env bash
# =============================================================================
# AbacusCopilot server upgrade helper
# =============================================================================
# Swaps the current install's source directory for a new release tarball in
# place, refreshes the editable pip registration, and preserves the platform
# Bader binary. Does NOT need internet or setup.sh: the conda env and all
# Python dependencies are reused as-is (setups are editable `pip install -e .`).
#
# Usage:
#   bash upgrade.sh /path/to/abacuscopilot_v0.1.34_20260902.tar.gz
#   bash upgrade.sh                          # picks the newest *.tar.gz in cwd
#   bash upgrade.sh --install-deps <tarball> # also pip-install deps (needs a
#                                            # reachable PyPI mirror); use when
#                                            # the release changed pyproject deps
#
# Behavior:
#   1. Verify the tarball, locate the current install dir from the LIVE
#      editable install in the conda env.
#   2. Back up the current dir to <dir>_old_bak_<timestamp> (keep until you've
#      verified the new version, then delete).
#   3. Extract the new tarball to the SAME path (so the editable link still
#      points at it).
#   4. Restore extra libraries the user had dropped into the old install's
#      PP-Orb/ — the release tarball ships only the bundled default series
#      (SG15 + lanthanides), so external series downloaded by hand (e.g.
#      Dojo-NC-FR, ABACUS-APNS-PPORBs-v1) are merged back from the backup.
#   5. Preserve scripts/bader/bader.x — the tarball deliberately excludes the
#      platform binary; it would otherwise be lost on a full dir swap.
#   6. Refresh pip metadata with an offline reinstall, verify the version.
#
# Notes:
#   - `~/.abacuscopilot/config.yaml` lives outside the source dir → survives.
#   - If the release changed Python dependencies, run `bash setup.sh` instead
#     (now offline-safe), or pass --install-deps with a reachable mirror.
#   - Put this script somewhere stable on the server (e.g. ~/softwares/upgrade.sh).
#     It also ships inside every release tarball under scripts/.
# =============================================================================

set -euo pipefail

ENV_NAME="${ABACUS_ENV:-abacuscopilot}"
INSTALL_DEPS=0
TARBALL=""

for arg in "$@"; do
    case "$arg" in
        --install-deps) INSTALL_DEPS=1 ;;
        -h|--help)
            sed -n '2,35p' "$0"
            exit 0
            ;;
        -*) ;;
        *) TARBALL="$arg" ;;
    esac
done

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

echo ""
echo -e "  ${CYAN}${BOLD}AbacusCopilot Upgrade${NC}"

# ---- 0. Locate the tarball -------------------------------------------------
if [ -z "$TARBALL" ]; then
    TARBALL="$(ls -t abacuscopilot_v*.tar.gz 2>/dev/null | head -1 || true)"
fi
if [ -z "$TARBALL" ] || [ ! -f "$TARBALL" ]; then
    echo -e "  ${RED}Error: no tarball given or found.${NC}"
    echo -e "  Usage: bash upgrade.sh <abacuscopilot_v*.tar.gz>"
    exit 1
fi
TARBALL="$(cd "$(dirname "$TARBALL")" && pwd)/$(basename "$TARBALL")"
echo -e "  Tarball : ${BOLD}$(basename "$TARBALL")${NC}"

# ---- 1. Verify tarball -----------------------------------------------------
if ! gzip -t "$TARBALL" 2>/dev/null; then
    echo -e "  ${RED}Error: $TARBALL is not a valid gzip archive.${NC}"
    exit 1
fi
# Note: `|| true` — under `set -o pipefail`, `head -1` closing the pipe early
# makes tar die with SIGPIPE (141) and would abort the script silently.
TOPDIR="$(tar tzf "$TARBALL" 2>/dev/null | head -1 | cut -d/ -f1 || true)"
[ -n "$TOPDIR" ] || { echo -e "  ${RED}Error: empty archive.${NC}"; exit 1; }

# ---- 2. conda env must exist ------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
    echo -e "  ${RED}Error: 'conda' not found in PATH.${NC}"
    exit 1
fi
if ! conda run -n "$ENV_NAME" python -c "import sys" >/dev/null 2>&1; then
    echo -e "  ${RED}Error: conda env '${ENV_NAME}' not found.${NC}"
    echo -e "  Run ./setup.sh once on this server first, or set ABACUS_ENV."
    exit 1
fi

# ---- 3. Detect current install dir from the live editable install -----------
# abacuscopilot/__init__.py lives at <install>/abacuscopilot/__init__.py
PKG_ROOT="$(conda run -n "$ENV_NAME" python -c \
  "import abacuscopilot, os; print(os.path.dirname(os.path.dirname(os.path.abspath(abacuscopilot.__file__))))" \
  2>/dev/null || true)"
if [ -z "$PKG_ROOT" ] || [ ! -d "$PKG_ROOT" ]; then
    echo -e "  ${YELLOW}Live install not locatable — defaulting to ./abacuscopilot${NC}"
    PKG_ROOT="$(pwd)/abacuscopilot"
fi
PKG_ROOT="$(cd "$PKG_ROOT" 2>/dev/null && pwd || echo "$PKG_ROOT")"
[ -d "$PKG_ROOT" ] || { echo -e "  ${RED}Error: install dir $PKG_ROOT not found.${NC}"; exit 1; }
PARENT="$(dirname "$PKG_ROOT")"
echo -e "  Install : ${CYAN}${PKG_ROOT}${NC}"
echo -e "  Conda env: ${CYAN}${ENV_NAME}${NC}"

# ---- 4. Backup current install ----------------------------------------------
TS="$(date +%Y%m%d_%H%M%S)"
BACKUP="${PKG_ROOT}_old_bak_${TS}"
echo -e "  Backup  -> ${BOLD}${BACKUP}${NC}"
mv "$PKG_ROOT" "$BACKUP"

# ---- 5. Extract new tarball to the SAME location ----------------------------
echo -e "  Extracting new release..."
cd "$PARENT"
tar xzf "$TARBALL"
# Recreate the original path even if the tarball's top dir name differs.
if [ "$PARENT/$TOPDIR" != "$PKG_ROOT" ]; then
    mv "$PARENT/$TOPDIR" "$PKG_ROOT"
fi

# ---- 6. Restore extra libraries the user dropped into the old PP-Orb/ --------
# The release tarball deliberately ships only the bundled default series
# (SG15 + lanthanides + PP-Orb/README.md). Any external family the user
# downloaded into the previous install's PP-Orb/ (e.g. Dojo-NC-FR for the
# library-family auto-detection, ABACUS-APNS-PPORBs-v1, custom dirs) is NOT in
# the new release — merge it back from the backup. Only top-level entries that
# are MISSING from the new PP-Orb/ are copied, so the bundled SG15 / lanthanide
# trees that ship in the tarball stay authoritative (never overwritten or
# duplicated by an older backup copy).
if [ -d "$BACKUP/PP-Orb" ]; then
    RESTORED=0
    for entry in "$BACKUP"/PP-Orb/*; do
        [ -e "$entry" ] || continue
        name="$(basename "$entry")"
        if [ ! -e "$PKG_ROOT/PP-Orb/$name" ]; then
            mkdir -p "$PKG_ROOT/PP-Orb"
            cp -a "$entry" "$PKG_ROOT/PP-Orb/$name"
            RESTORED=$((RESTORED + 1))
        fi
    done
    if [ "$RESTORED" -gt 0 ]; then
        echo -e "        ${GREEN}✓${NC} restored ${RESTORED} extra PP-Orb library/libraries from the old install"
    fi
fi

# ---- 7. Preserve platform-compiled Bader binary ------------------------------
if [ -x "$BACKUP/scripts/bader/bader.x" ]; then
    mkdir -p "$PKG_ROOT/scripts/bader"
    cp -f "$BACKUP/scripts/bader/bader.x" "$PKG_ROOT/scripts/bader/bader.x"
    chmod +x "$PKG_ROOT/scripts/bader/bader.x"
    echo -e "        ${GREEN}✓${NC} preserved platform bader binary"
fi

# ---- 8. Refresh editable pip registration ------------------------------------
cd "$PKG_ROOT"
if [ "$INSTALL_DEPS" = "1" ]; then
    echo -e "  Installing with dependencies (needs a reachable PyPI mirror)..."
    conda run -n "$ENV_NAME" pip install -e . --no-build-isolation --upgrade
else
    echo -e "  Refreshing editable install (no deps, no network needed)..."
    conda run -n "$ENV_NAME" pip install -e . --no-build-isolation --no-deps --upgrade
fi

# ---- 9. Verify ---------------------------------------------------------------
echo -e "  Verifying..."
VER="$(cd "$PKG_ROOT" && conda run -n "$ENV_NAME" python -c \
  "import abacuscopilot; print('v' + abacuscopilot.__version__)" 2>&1 || echo "unknown")"
echo -e "        ${GREEN}✓${NC} Installed version: ${BOLD}${VER}${NC}"

EXPECTED="$(basename "$TARBALL" | sed -E 's/^abacuscopilot_(v[0-9][0-9.]*)_.*/\1/')"
if [ -n "$EXPECTED" ] && [ "$EXPECTED" != "$VER" ]; then
    echo -e "  ${YELLOW}! Tarball name says ${EXPECTED}, but installed version is ${VER}.${NC}"
    echo -e "  ${YELLOW}  Did you bump __init__.py / pyproject.toml before packaging?${NC}"
fi

echo ""
echo -e "  ${GREEN}${BOLD}Upgrade complete!${NC}"
echo -e "  Start:  conda activate ${ENV_NAME} && abacuscopilot"
echo -e "  Old install kept at: ${BOLD}${BACKUP}${NC}  (delete after verifying)"
echo -e "  ${YELLOW}Config ~/.abacuscopilot/config.yaml was NOT touched.${NC}"
echo ""
