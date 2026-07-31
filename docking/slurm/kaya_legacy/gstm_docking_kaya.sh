#!/bin/bash --login

#############################################
# SLURM directives
#SBATCH --job-name=gstm_docking
#SBATCH --partition=work
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err   

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOBID"
echo "Node: $SLURM_NODELIST"

#############################################
# Load required modules 
module load Anaconda3/2024.06
conda activate /group/sms056/conda/envs/vina

#############################################
# Define paths
EXECUTABLE=$SLURM_SUBMIT_DIR/Vina_rigid_M.pl
SCRATCH=$MYSCRATCH/$SLURM_JOB_NAME/$SLURM_JOBID
RESULTS=$SLURM_SUBMIT_DIR/results/$SLURM_JOBID

#############################################
# Validate inputs
if [ ! -f "$EXECUTABLE" ]; then
    echo "ERROR: $EXECUTABLE not found"
    exit 1
fi
if [ ! -f "$SLURM_SUBMIT_DIR/config_M1.txt" ]; then
    echo "ERROR: config_M1.txt not found"
    exit 1
fi
if [ ! -f "$SLURM_SUBMIT_DIR/ligand.txt" ]; then
    echo "ERROR: ligand.txt not found"
    exit 1
fi

#############################################
# Set up scratch
mkdir -p $SCRATCH $RESULTS

cp "$EXECUTABLE" "$SCRATCH/"
cp "$SLURM_SUBMIT_DIR/config_M1.txt" "$SCRATCH/"
cp "$SLURM_SUBMIT_DIR/ligand.txt" "$SCRATCH/"
cp "$SLURM_SUBMIT_DIR"/*.pdbqt "$SCRATCH/"

cd $SCRATCH

#############################################
# Run
printf "ligand.txt\n" | perl Vina_rigid_M.pl
EXIT_CODE=$?

echo "Finished at: $(date) — exit code: $EXIT_CODE"

#############################################
# Copy results to permanent storage
cp -r $SCRATCH/* $RESULTS/

#############################################
# Clean up scratch
rm -rf $SCRATCH
cd $HOME

echo "Results saved to: $RESULTS"
[ $EXIT_CODE -ne 0 ] && exit $EXIT_CODE
