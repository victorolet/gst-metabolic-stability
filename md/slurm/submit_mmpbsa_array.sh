#!/bin/bash
# Submission helper for run_mmpbsa_array.sh -- computes the --array=X-Y
# range (has to be known at sbatch time, can't be derived inside the batch
# script) from the same frame args and submits for you.
#
# Plain bash -- run directly, do NOT sbatch it.
#
# Usage (from pipeline_development/):
#   bash md/slurm/submit_mmpbsa_array.sh <run-dir> <start-frame> <end-frame> [interval] [frames-per-chunk] [max-concurrent]
# Example (140 sampled frames, chunks of 20 -> 7 array tasks, max 8 running at once):
#   bash md/slurm/submit_mmpbsa_array.sh md/runs/gsta1_lig03506 300 1000 5 20 8
#
# After all tasks finish:
#   python3 md/tools/merge_mmpbsa_chunks.py md/runs/gsta1_lig03506/mmpbsa_chunks
set -euo pipefail

RUN_DIR=${1:-}
STARTFRAME=${2:-}
ENDFRAME=${3:-}
INTERVAL=${4:-1}
CHUNK_SIZE=${5:-20}
MAX_CONCURRENT=${6:-8}

if [ -z "$RUN_DIR" ] || [ -z "$STARTFRAME" ] || [ -z "$ENDFRAME" ]; then
    echo "ERROR: usage: bash submit_mmpbsa_array.sh <run-dir> <start-frame> <end-frame> [interval] [frames-per-chunk] [max-concurrent]"
    exit 1
fi

TOTAL_FRAMES=$(( (ENDFRAME - STARTFRAME) / INTERVAL ))
if [ "$TOTAL_FRAMES" -le 0 ]; then
    echo "ERROR: computed 0 sampled frames from start=$STARTFRAME end=$ENDFRAME interval=$INTERVAL"
    exit 1
fi
N_CHUNKS=$(( (TOTAL_FRAMES + CHUNK_SIZE - 1) / CHUNK_SIZE ))

echo "Total sampled frames: $TOTAL_FRAMES   Chunk size: $CHUNK_SIZE   -> $N_CHUNKS array tasks (0-$((N_CHUNKS - 1))), max $MAX_CONCURRENT concurrent"

mkdir -p md/logs

sbatch --array=0-$((N_CHUNKS - 1))%${MAX_CONCURRENT} \
    md/slurm/run_mmpbsa_array.sh "$RUN_DIR" "$STARTFRAME" "$ENDFRAME" "$INTERVAL" "$CHUNK_SIZE"
