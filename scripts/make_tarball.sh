#!/usr/bin/env bash
# =============================================================================
# AbacusCopilot release tarball builder
# =============================================================================
# Builds abacuscopilot_v{version}_{date}.tar.gz at the repo's parent dir,
# bundling the source + the PP-Orb library tree (for out-of-the-box use).
#
# PP-Orb/ ships the bundled default SG15 (+ lanthanide) trees; the
# PP-Orb/README.md placeholder (tracked in git) stays in the archive to guide
# users who add their own pseudopotential/orbital libraries.
#
# Any external series the user registered under libraries.families (APNS,
# Dojo-NC-FR, ...) and downloaded into PP-Orb/ is NOT bundled into the release
# (see the "External library families" README section) — those folders are
# excluded below so they never inflate or leak into the tarball.
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
  # External user-registered series downloaded under PP-Orb/ are not bundled
  # (kept out of the release tarball — the config libraries.families entries
  # point at them locally; see README "External library families").
  --exclude='PP-Orb/ABACUS-APNS-PPORBs-v1'
  --exclude='PP-Orb/Dojo-NC-FR'
  # Stray archive bundles left at the repo root (e.g. a developer-zipped
  # PP-Orb.zip / PP-Orb.7z for a manual server transfer) must not ride along
  # in the release.
  --exclude='*.zip'
  --exclude='*.7z'
  --exclude='*.tar.gz'
  --exclude='*.tgz'
)

rm -f "$OUT"
# macOS bsdtar stores Apple extended attributes (com.apple.quarantine on
# downloaded files etc.) as PAX headers that GNU tar on Linux servers prints
# hundreds of "Ignoring unknown extended header keyword" warnings for on
# extract. Both flags together are needed to strip them from the archive.
tar czf "$OUT" --no-xattrs --no-mac-metadata "${EXCLUDES[@]}" "$(basename "$ROOT")"
gzip -t "$OUT"

echo "OK: $(ls -lh "$OUT" | awk '{print $5}') tarball written"
