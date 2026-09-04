#!/usr/bin/env bash
# =============================================================================
# AbacusCopilot release tarball builder
# =============================================================================
# Builds abacuscopilot_v{version}_{date}.tar.gz at the repo's parent dir,
# bundling the source + the PP-Orb library tree (for out-of-the-box use).
#
# Excluded from the tarball (users don't need these):
#   - .git, editor/IDE + OS caches, Claude Code dirs, egg-info
#   - test/ (local SLURM/DeepMD job scripts — machine-specific)
#   - tests/ (pytest suite — for developers only; still in git)
#   - compiled platform-specific bader binary (scripts/bader/bader.x;
#     setup.sh builds/downloads it on the target machine)
#
# Usage:
#   bash scripts/make_tarball.sh            # builds abacuscopilot_vX.Y.Z_YYYYMMDD.tar.gz
# =============================================================================

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # repo root
PARENT="$(dirname "$ROOT")"                              # dir that holds the repo dir
cd "$PARENT"

VERSION=$(grep -m1 '__version__ *=' "$ROOT/abacuscopilot/__init__.py" | sed -E 's/.*"([^"]+)".*/\1/')
VDATE=$(grep -m1 '__version_date__ *=' "$ROOT/abacuscopilot/__init__.py" | sed -E 's/.*"([^"]+)".*/\1/')
OUT="$PARENT/abacuscopilot_v${VERSION}_${VDATE//-/}.tar.gz"

echo "Packaging AbacusCopilot v${VERSION} (${VDATE}) -> ${OUT}"

EXCLUDES=(
  --exclude='.git'
  --exclude='.claude'
  --exclude='.agents'
  --exclude='server-profiles'
  --exclude='*.egg-info'
  --exclude='.coverage'
  --exclude='.pytest_cache'
  --exclude='.mypy_cache'
  --exclude='.ruff_cache'
  --exclude='__pycache__'
  --exclude='.DS_Store'
  --exclude='test'
  --exclude='tests'
  --exclude='scripts/bader/bader.x'
)

rm -f "$OUT"
# macOS bsdtar stores Apple extended attributes (com.apple.quarantine on
# downloaded files etc.) as PAX headers that GNU tar on Linux servers prints
# hundreds of "Ignoring unknown extended header keyword" warnings for on
# extract. Both flags together are needed to strip them from the archive.
tar czf "$OUT" --no-xattrs --no-mac-metadata "${EXCLUDES[@]}" "$(basename "$ROOT")"
gzip -t "$OUT"

echo "OK: $(ls -lh "$OUT" | awk '{print $5}') tarball written"
