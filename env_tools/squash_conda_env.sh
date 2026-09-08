#!/bin/bash
# DEPRECATED on Setonix -- confirmed 2026-09-07: mksquashfs builds fine,
# but the mount step (squashfuse, unprivileged FUSE) fails with "Operation
# not permitted". Setonix blocks unprivileged FUSE mounts at the platform
# level (Apptainer's own mounting works because it uses a privileged,
# admin-installed helper a personal spack-built squashfuse doesn't get).
# Kept for reference only -- use env_tools/backup_env_to_acacia.sh +
# scratch instead (see docs/implementation_notes.md).
#
# Pack a conda environment into a single SquashFS image, to get it out of
# Pawsey's /software file-count quota (a mature env's tens of thousands of
# files become ONE file on disk). Smoke-tests the image via a temporary
# mount but does NOT touch the original environment -- deleting/swapping it
# in is a separate manual step, printed at the end, so you can verify first.
#
# Run this ON Setonix (needs mksquashfs + squashfuse). See
# ../docs/implementation_notes.md for why this exists and how mounting
# works day-to-day (see mount_conda_env.sh).
#
# Usage:
#   bash squash_conda_env.sh <env_dir> <out_sqfs_path>
# Example:
#   bash squash_conda_env.sh /software/projects/pawsey1376/volet/conda_envs/gst_ml \
#       /software/projects/pawsey1376/volet/squashed_envs/gst_ml.sqfs
set -euo pipefail

ENV_DIR=${1:-}
OUT_SQFS=${2:-}

if [ -z "$ENV_DIR" ] || [ -z "$OUT_SQFS" ]; then
    echo "Usage: $0 <env_dir> <out_sqfs_path>"
    exit 1
fi
if [ ! -d "$ENV_DIR" ]; then
    echo "ERROR: $ENV_DIR not found"
    exit 1
fi
if [ ! -x "$ENV_DIR/bin/python3" ] && [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "WARNING: no bin/python(3) found under $ENV_DIR -- is this really a conda env directory?"
fi

command -v mksquashfs >/dev/null 2>&1 || { echo "ERROR: mksquashfs not found. Try: module load squashfs (confirmed on Setonix as squashfs/4.6.1)"; exit 1; }
command -v squashfuse >/dev/null 2>&1 || { echo "ERROR: squashfuse not found. Try: module load squashfuse, or module spider squashfuse (ask Pawsey support if no such module)"; exit 1; }

ENV_DIR=$(cd "$ENV_DIR" && pwd)
ENV_NAME=$(basename "$ENV_DIR")
mkdir -p "$(dirname "$OUT_SQFS")"

FILE_COUNT=$(find "$ENV_DIR" | wc -l)
echo "======================================================================"
echo "STAGE 1: building squashfs image"
echo "======================================================================"
echo "$ENV_DIR currently contains $FILE_COUNT files/dirs."
mksquashfs "$ENV_DIR" "$OUT_SQFS" -noappend -processors 8

echo ""
echo "======================================================================"
echo "STAGE 2: smoke-test via a temporary mount"
echo "======================================================================"
TEST_MOUNT=$(mktemp -d)
cleanup() { fusermount -u "$TEST_MOUNT" 2>/dev/null; rmdir "$TEST_MOUNT" 2>/dev/null || true; }
trap cleanup EXIT

squashfuse "$OUT_SQFS" "$TEST_MOUNT"

PYBIN="$TEST_MOUNT/bin/python3"
[ -x "$PYBIN" ] || PYBIN="$TEST_MOUNT/bin/python"
if [ -x "$PYBIN" ]; then
    echo "[python] $("$PYBIN" --version 2>&1)"
else
    echo "WARNING: no python binary found in the mounted image -- inspect $TEST_MOUNT manually."
fi
echo "[mount check] $(ls "$TEST_MOUNT" | wc -l) top-level entries visible in the mounted image"

echo ""
echo "======================================================================"
echo "[OK] $OUT_SQFS built and mounts cleanly ($(du -h "$OUT_SQFS" | cut -f1))"
echo "======================================================================"
echo ""
echo "This script has NOT touched $ENV_DIR. Once you've verified the mounted"
echo "image actually works for real work (activate it from $TEST_MOUNT-style"
echo "path and run something from $ENV_NAME before it's unmounted, or repeat"
echo "this smoke test), swap it in manually:"
echo ""
echo "  1. Make sure nothing has $ENV_NAME active right now."
echo "  2. rm -rf \"$ENV_DIR\""
echo "  3. mkdir -p \"$ENV_DIR\""
echo "  4. source env_tools/mount_conda_env.sh && mount_conda_env \"$ENV_DIR\" \"$OUT_SQFS\""
echo "  5. conda activate \"$ENV_DIR\"   # should now be backed by the squashfs image"
echo ""
echo "Reminder: squashfuse mounts don't survive a new login shell or a new"
echo "SLURM job's node allocation -- mount_conda_env.sh handles re-mounting"
echo "automatically and is already wired into the SLURM scripts that use"
echo "gst_ml / gmxMMPBSA, so this only matters if you activate interactively."
