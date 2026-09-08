#!/bin/bash --login
# Batch wrapper around ../03_run_md.sh (production MD is hours long, needs
# to run as a batch job not an interactive session). Validated performance:
# ~43 ns/day (1 node, 32 OpenMP threads, CPU-only). See
# docs/implementation_notes.md for background.
#
# ONE-TIME SETUP: mkdir -p md/logs
#
# Usage (submit from pipeline_development/):
#   sbatch md/slurm/run_md.sh <protein-dir> <gsh-dir> <ligand-dir> <outdir>
# Example:
#   sbatch md/slurm/run_md.sh \
#       md/system/gsta1_prepped md/system/gsta1_gsh_prepped \
#       md/system/gsta1_lig03506_prepped md/runs/gsta1_lig03506

#SBATCH --account=pawsey1376
#SBATCH --job-name=gst_md
#SBATCH --partition=work
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=md/logs/%x_%j.out
#SBATCH --error=md/logs/%x_%j.err

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID  Node: $SLURM_NODELIST"

#############################################
# Args
PROTEIN_DIR=$1
GSH_DIR=$2
LIGAND_DIR=$3
OUTDIR=$4
MDP_DIR=$SLURM_SUBMIT_DIR/md/mdp
FFDIR=$SLURM_SUBMIT_DIR/topology/charmm36-feb2026_cgenff-5.0.ff

if [ -z "$PROTEIN_DIR" ] || [ -z "$GSH_DIR" ] || [ -z "$LIGAND_DIR" ] || [ -z "$OUTDIR" ]; then
    echo "ERROR: usage: sbatch run_md.sh <protein-dir> <gsh-dir> <ligand-dir> <outdir>"
    exit 1
fi

#############################################
# Environment. gst_ml lives on scratch (moved off /software's file-count
# quota -- SquashFS/FUSE doesn't work on Setonix, see
# docs/implementation_notes.md). Back it up periodically with
# env_tools/backup_env_to_acacia.sh as insurance against scratch's 21-day
# purge.
source $MYSOFTWARE/miniconda3/etc/profile.d/conda.sh
conda activate $MYSCRATCH/conda_envs/gst_ml
module load gromacs/2024.3-mixed

#############################################
# Clean any partial output from a prior interactive/Ctrl-C'd attempt.
rm -rf "$OUTDIR"

bash "$SLURM_SUBMIT_DIR/md/03_run_md.sh" \
    --protein-dir "$PROTEIN_DIR" \
    --gsh-dir "$GSH_DIR" \
    --ligand-dir "$LIGAND_DIR" \
    --mdp-dir "$MDP_DIR" \
    --ffdir "$FFDIR" \
    --outdir "$OUTDIR"
EXIT_CODE=$?

echo "Finished at: $(date) -- exit code: $EXIT_CODE"
exit $EXIT_CODE
