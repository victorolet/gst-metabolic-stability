#!/bin/bash
# DEPRECATED on Setonix -- see squash_conda_env.sh's header (unprivileged
# FUSE mounts are blocked platform-wide; "Operation not permitted"
# confirmed 2026-09-07). No longer wired into any SLURM script. Kept for
# reference only -- use scratch + backup_env_to_acacia.sh instead.
#
# Idempotent mount-before-activate helper for a SquashFS-packed conda env.
# Source this (not `bash`-execute it) and call mount_conda_env, then
# `conda activate` as normal -- see ../docs/implementation_notes.md.
#
# Safe to add to a script unconditionally, whether or not an env has
# actually been converted yet: it only mounts when the .sqfs image exists
# AND env_dir is currently empty/missing (a live directory with real files
# is never shadowed) -- so it acts as a fallback/restore (e.g. after a
# scratch purge wipes a live env) rather than something that could hide a
# live copy you're still using directly.
#
# Usage (source it, then call the function):
#   source env_tools/mount_conda_env.sh
#   mount_conda_env /software/projects/pawsey1376/volet/conda_envs/gst_ml \
#       /software/projects/pawsey1376/volet/squashed_envs/gst_ml.sqfs
#   conda activate /software/projects/pawsey1376/volet/conda_envs/gst_ml

mount_conda_env() {
    local env_dir=$1
    local sqfs=$2

    if [ -z "$env_dir" ] || [ -z "$sqfs" ]; then
        echo "mount_conda_env: usage: mount_conda_env <env_dir> <sqfs_path>" >&2
        return 1
    fi
    if [ ! -f "$sqfs" ]; then
        # Not squashed yet -- nothing to do, env_dir is used live.
        return 0
    fi
    if mountpoint -q "$env_dir" 2>/dev/null; then
        # Already mounted (e.g. function called twice in the same session).
        return 0
    fi
    if [ -d "$env_dir" ] && [ -n "$(ls -A "$env_dir" 2>/dev/null)" ]; then
        # Live directory still has real content -- use it as-is, don't shadow it.
        return 0
    fi

    mkdir -p "$env_dir"
    if ! command -v squashfuse >/dev/null 2>&1; then
        echo "mount_conda_env: ERROR: $sqfs exists but squashfuse is not on PATH (try: module load squashfuse)" >&2
        return 1
    fi
    squashfuse "$sqfs" "$env_dir"
    echo "mount_conda_env: mounted $sqfs -> $env_dir"
}
