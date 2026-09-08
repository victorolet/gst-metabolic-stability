#!/bin/bash
# Restore a conda env directory from an Acacia backup -- e.g. after
# scratch's 21-day purge wipes a live copy. Always restore to the SAME
# path the env was originally built at: conda environments bake absolute
# paths into scripts/shebangs/conda-meta, so a restore to a different path
# can silently break (see docs/implementation_notes.md). Refuses to
# overwrite a directory that already has real content.
#
# Usage:
#   bash restore_env_from_acacia.sh <remote>:<bucket>/<path> <env_dir>
# Example:
#   bash restore_env_from_acacia.sh acacia:volet-envs/gst_ml $MYSCRATCH/conda_envs/gst_ml
set -euo pipefail

REMOTE_PATH=${1:-}
ENV_DIR=${2:-}

if [ -z "$REMOTE_PATH" ] || [ -z "$ENV_DIR" ]; then
    echo "Usage: $0 <remote>:<bucket>/<path> <env_dir>"
    exit 1
fi
command -v rclone >/dev/null 2>&1 || { echo "ERROR: rclone not found. Try: module load rclone"; exit 1; }

if [ -d "$ENV_DIR" ] && [ -n "$(ls -A "$ENV_DIR" 2>/dev/null)" ]; then
    echo "ERROR: $ENV_DIR already has content -- refusing to overwrite."
    echo "Move or remove it first if you really mean to restore over it."
    exit 1
fi

mkdir -p "$ENV_DIR"
echo "Restoring $REMOTE_PATH -> $ENV_DIR"
rclone sync "$REMOTE_PATH" "$ENV_DIR" --progress

echo ""
echo "[OK] Restored. Verify before trusting it for real work:"
echo "  source \$MYSOFTWARE/miniconda3/etc/profile.d/conda.sh"
echo "  conda activate \"$ENV_DIR\""
echo "  python3 -c \"import sys; print(sys.version)\""
