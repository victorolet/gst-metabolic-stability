#!/bin/bash --login
# MM-PBSA rescoring -- SLURM ARRAY version. Real MPI parallelism is blocked
# on Setonix (GROMACS/mpi4py library mismatch -- see
# docs/implementation_notes.md), so this splits the frame range into
# independent chunks and runs one fully serial gmx_MMPBSA per chunk
# (--nranks 1, the validated config) as array tasks -- embarrassingly
# parallel, no MPI risk.
#
# Do not sbatch this directly with a guessed --array range -- use
# submit_mmpbsa_array.sh, which computes it for you.
#
# Each task writes to $RUN_DIR/mmpbsa_chunks/chunk_<NN>/. Once all tasks
# finish, merge with:
#   python3 md/tools/merge_mmpbsa_chunks.py <run-dir>/mmpbsa_chunks

#SBATCH --account=pawsey1376
#SBATCH --job-name=gst_mmpbsa_arr
#SBATCH --partition=work
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=md/logs/%x_%A_%a.out
#SBATCH --error=md/logs/%x_%A_%a.err

echo "Array task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID started at: $(date)"

#############################################
# Args (identical for every array task -- each task derives its own slice below)
RUN_DIR=$1
G_STARTFRAME=$2
G_ENDFRAME=$3
INTERVAL=${4:-1}
CHUNK_SIZE=${5:-20}   # sampled frames per chunk, matches the validated serial batch

if [ -z "$RUN_DIR" ] || [ -z "$G_STARTFRAME" ] || [ -z "$G_ENDFRAME" ]; then
    echo "ERROR: usage: sbatch --array=0-N-1 run_mmpbsa_array.sh <run-dir> <start-frame> <end-frame> [interval] [frames-per-chunk]"
    echo "        (normally you don't call this directly -- use submit_mmpbsa_array.sh)"
    exit 1
fi
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "ERROR: this script must be submitted with --array=0-N-1 set. Use submit_mmpbsa_array.sh instead of calling sbatch on this directly."
    exit 1
fi

#############################################
# Derive this task's own [chunk-start, chunk-end) frame slice.
CHUNK_SPAN=$(( CHUNK_SIZE * INTERVAL ))   # span in original-trajectory frame units
STARTFRAME=$(( G_STARTFRAME + SLURM_ARRAY_TASK_ID * CHUNK_SPAN ))
ENDFRAME=$(( STARTFRAME + CHUNK_SPAN ))
if [ "$ENDFRAME" -gt "$G_ENDFRAME" ]; then
    ENDFRAME=$G_ENDFRAME
fi
if [ "$STARTFRAME" -ge "$G_ENDFRAME" ]; then
    echo "Task $SLURM_ARRAY_TASK_ID: chunk start ($STARTFRAME) is past the global end frame ($G_ENDFRAME) -- nothing to do, exiting cleanly."
    exit 0
fi

CHUNK_OUT="$RUN_DIR/mmpbsa_chunks/chunk_$(printf '%02d' $SLURM_ARRAY_TASK_ID)"
echo "Task $SLURM_ARRAY_TASK_ID: frames [$STARTFRAME, $ENDFRAME) interval $INTERVAL -> $CHUNK_OUT"

#############################################
# Environment -- see run_mmpbsa.sh for the scratch + Acacia-backup pattern.
source $MYSOFTWARE/miniconda3/etc/profile.d/conda.sh
conda activate $MYSCRATCH/conda_envs/gmxMMPBSA
module load gromacs/2024.3-mixed

#############################################
# Run (serial -- --nranks defaults to 1 in 05_run_mmpbsa.sh, the validated config)
bash "$SLURM_SUBMIT_DIR/md/05_run_mmpbsa.sh" \
    --run-dir "$RUN_DIR" \
    --out-dir "$CHUNK_OUT" \
    --start-frame "$STARTFRAME" \
    --end-frame "$ENDFRAME" \
    --interval "$INTERVAL"
EXIT_CODE=$?

echo "Task $SLURM_ARRAY_TASK_ID finished at: $(date) -- exit code: $EXIT_CODE"
exit $EXIT_CODE
