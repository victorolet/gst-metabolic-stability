#!/bin/bash
# Back up a conda env directory to Acacia (Pawsey's S3-compatible object
# storage) via rclone. Acacia has a space quota, not a file-count quota, so
# it's fine for a many-small-files conda env -- unlike /software. Use this
# periodically for envs living on scratch (21-day purge) as insurance; see
# docs/implementation_notes.md.
#
# One-time setup:
#   1. Generate an Acacia access key via the Pawsey portal.
#   2. Add a stanza to ~/.config/rclone/rclone.conf -- see
#      acacia_rclone.conf.example in this directory.
#   3. Create the bucket once: rclone mkdir acacia:<bucket>
#
# Uses --copy-links: conda envs are full of symlinks (bin/python -> a
# versioned binary, shared-library links, etc.), and rclone skips symlinks
# by default -- --copy-links dereferences them and backs up the real
# content instead, so the backup is actually complete and self-contained
# (a plain sync would silently produce a backup missing every symlinked
# file). Costs some extra space vs a true symlink-preserving backup.
#
# Defaults to a dry run (prints what WOULD change, syncs nothing). Pass
# --do-it as the third argument to actually upload.
#
# Usage:
#   bash backup_env_to_acacia.sh <env_dir> <remote>:<bucket>/<path> [--do-it]
# Example:
#   bash backup_env_to_acacia.sh $MYSCRATCH/conda_envs/gst_ml acacia:volet-envs/gst_ml --do-it
set -euo pipefail

ENV_DIR=${1:-}
REMOTE_PATH=${2:-}
DO_IT=${3:-}

if [ -z "$ENV_DIR" ] || [ -z "$REMOTE_PATH" ]; then
    echo "Usage: $0 <env_dir> <remote>:<bucket>/<path> [--do-it]"
    exit 1
fi
if [ ! -d "$ENV_DIR" ]; then
    echo "ERROR: $ENV_DIR not found"
    exit 1
fi
command -v rclone >/dev/null 2>&1 || { echo "ERROR: rclone not found. Try: module load rclone"; exit 1; }

if [ "$DO_IT" = "--do-it" ]; then
    echo "Backing up $ENV_DIR -> $REMOTE_PATH"
    rclone sync "$ENV_DIR" "$REMOTE_PATH" --copy-links --progress
    echo "[OK] Backed up to $REMOTE_PATH"
else
    echo "DRY RUN (no files uploaded) -- $ENV_DIR -> $REMOTE_PATH"
    rclone sync "$ENV_DIR" "$REMOTE_PATH" --copy-links --dry-run
    echo ""
    echo "Review the above, then re-run with --do-it to actually upload:"
    echo "  bash $0 \"$ENV_DIR\" \"$REMOTE_PATH\" --do-it"
fi
