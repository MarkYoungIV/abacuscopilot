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
#      A top-level library the new layout has since folded into a series folder
#      (the old flat `SG15-Version1p0_Pseudopotential` and friends) is left out
#      — the series now ships it, and restoring the stray copy would put the
#      same series in the tree twice.
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
# A non-interactive shell does not source ~/.bashrc, so `ssh host 'bash
# upgrade.sh ...'` starts with no conda on PATH even on machines where it works
# fine interactively.  Look in the usual install locations before giving up.
if ! command -v conda >/dev/null 2>&1; then
    for _conda_sh in \
        "$HOME/miniconda3/etc/profile.d/conda.sh" \
        "$HOME/anaconda3/etc/profile.d/conda.sh" \
        "$HOME/miniconda/etc/profile.d/conda.sh" \
        /opt/conda/etc/profile.d/conda.sh
    do
        if [ -r "$_conda_sh" ]; then
            # shellcheck disable=SC1090
            . "$_conda_sh"
            break
        fi
    done
fi

if ! command -v conda >/dev/null 2>&1; then
    echo -e "  ${RED}Error: 'conda' not found in PATH.${NC}"
    echo -e "  Source it first, e.g. ${BOLD}. ~/miniconda3/etc/profile.d/conda.sh${NC}"
    exit 1
fi
if ! conda run -n "$ENV_NAME" python -c "import sys" >/dev/null 2>&1; then
    echo -e "  ${RED}Error: conda env '${ENV_NAME}' not found.${NC}"
    echo -e "  Run ./setup.sh once on this server first, or set ABACUS_ENV."
    exit 1
fi

# ---- 3. Detect current install dir from the live editable install -----------
# The install dir is the release root: the tree holding pyproject.toml, PP-Orb/
# and scripts/.  It is found by walking up from the package to that marker, not
# by counting dirname() calls — the package sits at <root>/src/abacuscopilot
# under the src layout and at <root>/abacuscopilot before it, so a fixed count
# silently returns <root>/src, and the script would then swap out the wrong
# tree.  The marker is layout-independent.
#
# `python -I` keeps the cwd off sys.path so the answer cannot depend on where
# the script was launched from.  Under the old flat layout the install dir's own
# parent held a directory named abacuscopilot, and importing from there resolved
# to it as a namespace package — __file__ came back None and this probe died in
# abspath.  The src layout removes the collision, but the probe should not rely
# on that.
PKG_ROOT="$(conda run -n "$ENV_NAME" python -I -c '
import os, abacuscopilot
p = os.path.dirname(os.path.abspath(abacuscopilot.__file__))
while not os.path.isfile(os.path.join(p, "pyproject.toml")):
    parent = os.path.dirname(p)
    if parent == p:
        p = ""
        break
    p = parent
print(p)
' 2>/dev/null | tail -1 || true)"
if [ -z "$PKG_ROOT" ] || [ ! -d "$PKG_ROOT" ]; then
    # Last resort, and only when the guess actually looks like an install — a
    # silent fallback to a directory that happens to exist is how a wrong tree
    # gets swapped out.
    if [ -f "$(pwd)/abacuscopilot/pyproject.toml" ]; then
        PKG_ROOT="$(pwd)/abacuscopilot"
        echo -e "  ${YELLOW}Probe failed; falling back to ./abacuscopilot — check the path below.${NC}"
    else
        echo -e "  ${RED}Error: cannot locate the live install directory.${NC}"
        echo -e "  Run this from the install dir's parent, or set ABACUS_ENV."
        exit 1
    fi
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
        # A top-level library that the NEW layout folded into a series folder is
        # not an extra library — it is that series' own old copy.  Without this,
        # an upgrade from the flat layout restores
        # `SG15-Version1p0_Pseudopotential` beside
        # `SG15-Version1p0/SG15-Version1p0_Pseudopotential`: the install then
        # holds one series twice, the stray copy belonging to none.  Skipped only
        # when the new location EXISTS, so a release that one day drops a series
        # cannot take away the user's only copy of it.  `case`, not an
        # associative array — macOS ships bash 3.2.
        case "$name" in
            SG15-Version1p0_Pseudopotential)
                superseded="$PKG_ROOT/PP-Orb/SG15-Version1p0/SG15-Version1p0_Pseudopotential" ;;
            SG15-Version1p0__StandardOrbitals-Version2p0)
                superseded="$PKG_ROOT/PP-Orb/SG15-Version1p0/SG15-Version1p0__StandardOrbitals-Version2p0" ;;
            PD04.3+f--core.icmod1)
                superseded="$PKG_ROOT/PP-Orb/lanthanides-f--core.icmod1/PD04.3+f--core.icmod1" ;;
            *)
                superseded="" ;;
        esac
        if [ -n "$superseded" ] && [ -e "$superseded" ]; then
            echo -e "        ${YELLOW}-${NC} ${name}: SKIPPED (superseded, now at PP-Orb/${superseded#"$PKG_ROOT/PP-Orb/"})"
            continue
        fi
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
# The version is read from the extracted tree, never parsed out of the tarball
# name. A name-parsing regex over [0-9.] looked fine but silently failed to match
# any letter-suffixed release (v0.1.35b, v0.1.35c): it left the comparison
# against the whole filename, so correctly packaged releases were reported as
# "did you forget to bump?" mismatches. The filename is only a label now — a
# stale one is worth a note, not an accusation.
echo -e "  Verifying..."
# -I again, so this reads the *registered* install rather than whatever the cwd
# happens to shadow it with — the point of the check is the editable link.
VER="$(cd "$PKG_ROOT" && conda run -n "$ENV_NAME" python -I -c \
  "import abacuscopilot; print(abacuscopilot.__version__)" 2>/dev/null | tail -1 || true)"
[ -n "$VER" ] || VER="unknown"
# What pip registered the install as — must agree with __init__.py.
PYVER="$(awk -F'"' '/^version[[:space:]]*=/{print $2; exit}' "$PKG_ROOT/pyproject.toml")"
[ -n "$PYVER" ] || PYVER="unknown"

echo -e "        ${GREEN}✓${NC} Installed version: ${BOLD}v${VER}${NC}"

if [ "$VER" != "$PYVER" ]; then
    echo -e "  ${YELLOW}! __init__.py says v${VER}, but pyproject.toml says v${PYVER}.${NC}"
    echo -e "  ${YELLOW}  Bump both before packaging.${NC}"
    exit 1
fi

# Advisory only: the archive's name is a convenience for the user, and a
# mismatch here says nothing about whether the upgrade itself worked.
NAMEVER="$(basename "$TARBALL" | sed -nE 's/^abacuscopilot_v?([0-9][0-9A-Za-z.]*)_.*$/\1/p')"
if [ -n "$NAMEVER" ] && [ "$NAMEVER" != "$VER" ]; then
    echo -e "  ${YELLOW}! Tarball name says v${NAMEVER}, the tree says v${VER}.${NC}"
    echo -e "  ${YELLOW}  Install is fine — only the archive's name is stale.${NC}"
fi

echo ""
echo -e "  ${GREEN}${BOLD}Upgrade complete!${NC}"
echo -e "  Start:  conda activate ${ENV_NAME} && abacuscopilot"
echo -e "  Old install kept at: ${BOLD}${BACKUP}${NC}  (delete after verifying)"
echo -e "  ${YELLOW}Config ~/.abacuscopilot/config.yaml was NOT touched.${NC}"
echo ""
