#!/bin/bash --login
# Batch wrapper around ../05_run_mmpbsa.sh with MPI parallel rescoring
# (--ntasks ranks, matched to --nranks so mpirun's slot count matches what
# SLURM granted). NOTE: on Setonix this MPI path is blocked by a
# GROMACS/mpi4py library mismatch -- see docs/implementation_notes.md and
# prefer run_mmpbsa_array.sh / submit_mmpbsa_array.sh instead. Kept for
# reference / in case that gets resolved.
#
# ONE-TIME SETUP: mkdir -p md/logs
#
# Usage (submit from pipeline_development/):
#   sbatch md/slurm/run_mmpbsa.sh <run-dir> <start-frame> <end-frame> [interval]

#SBATCH --account=pawsey1376
#SBATCH --job-name=gst_mmpbsa
#SBATCH --partition=work
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=md/logs/%x_%j.out
#SBATCH --error=md/logs/%x_%j.err

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID  Node: $SLURM_NODELIST  Ranks: $SLURM_NTASKS"

#############################################
# Args
RUN_DIR=$1
STARTFRAME=$2
ENDFRAME=$3
INTERVAL=${4:-1}

if [ -z "$RUN_DIR" ] || [ -z "$STARTFRAME" ] || [ -z "$ENDFRAME" ]; then
    echo "ERROR: usage: sbatch run_mmpbsa.sh <run-dir> <start-frame> <end-frame> [interval]"
    exit 1
fi

#############################################
# Environment -- gmxMMPBSA is its own conda env, kept separate from gst_ml.
# Lives on scratch (21-day purge); back it up periodically with
# env_tools/backup_env_to_acacia.sh as insurance. See
# docs/implementation_notes.md.
source $MYSOFTWARE/miniconda3/etc/profile.d/conda.sh
conda activate $MYSCRATCH/conda_envs/gmxMMPBSA
module load gromacs/2024.3-mixed

#############################################
# Run
bash "$SLURM_SUBMIT_DIR/md/05_run_mmpbsa.sh" \
    --run-dir "$RUN_DIR" \
    --start-frame "$STARTFRAME" \
    --end-frame "$ENDFRAME" \
    --interval "$INTERVAL" \
    --nranks "$SLURM_NTASKS"
EXIT_CODE=$?

echo "Finished at: $(date) -- exit code: $EXIT_CODE"
exit $EXIT_CODE
