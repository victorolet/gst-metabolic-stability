#!/bin/bash --login

#############################################
# AutoDock Vina docking -- GSTP1-1 -- Setonix (Pawsey)
#
# Same structure as gsta_docking.sh (see there for the full rationale),
# pointed at the GSTP1-1 receptor/config/perl script. Kaya version kept as
# gstp_docking_kaya.sh for reference.
#
# ONE-TIME SETUP before first use:
#   1. mkdir -p logs results
#   2. Create a `vina` conda env if one doesn't exist:
#        source $MYSOFTWARE/miniconda3/etc/profile.d/conda.sh
#        conda create -n vina -c conda-forge -c bioconda vina -y
#   3. Run prepare_ligands.py first to generate ligand_pdbqt/ and ligand.txt.
#
# CHUNK_SIZE / --time are first guesses -- submit a small array
# (--array=0-2) and check real per-chunk time via `sacct` before scaling up.
#
# Folder layout expected (relative to docking/):
#   receptors/1AQW_GSH.pdbqt  configs/config_P1.txt  scripts/Vina_rigid_P.pl
#   ligand_pdbqt/  ligand.txt  logs/  results/
#
# SUBMIT FROM THE docking/ DIRECTORY, array size derived from ligand count:
#   cd docking/
#   N=$(wc -l < ligand.txt); CH=50
#   sbatch --array=0-$(( (N + CH - 1) / CH - 1 )) slurm/gstp_docking.sh
#############################################

#SBATCH --account=pawsey1376
#SBATCH --job-name=gstp_docking
#SBATCH --partition=work
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

echo "Job started at: $(date)"
echo "Array job ID: $SLURM_ARRAY_JOB_ID  task: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"

#############################################
# Load Vina (Setonix / pawsey1376 conda pattern)
source $MYSOFTWARE/miniconda3/etc/profile.d/conda.sh
conda activate vina

#############################################
# Paths
EXECUTABLE=$SLURM_SUBMIT_DIR/scripts/Vina_rigid_P.pl
CONFIG=$SLURM_SUBMIT_DIR/configs/config_P1.txt
RECEPTOR=$SLURM_SUBMIT_DIR/receptors/1AQW_GSH.pdbqt
LIGAND_DIR=$SLURM_SUBMIT_DIR/ligand_pdbqt
MASTER_LIST=$SLURM_SUBMIT_DIR/ligand.txt
CHUNK_SIZE=50

SCRATCH=$MYSCRATCH/$SLURM_JOB_NAME/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
RESULTS=$SLURM_SUBMIT_DIR/results/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

#############################################
# Validate inputs
for f in "$EXECUTABLE" "$CONFIG" "$RECEPTOR" "$MASTER_LIST"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done
if [ ! -d "$LIGAND_DIR" ]; then
    echo "ERROR: $LIGAND_DIR not found -- run prepare_ligands.py first"
    exit 1
fi

#############################################
# Work out this task's slice of ligand.txt
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE + 1 ))
END=$(( START + CHUNK_SIZE - 1 ))
mkdir -p "$SCRATCH" "$RESULTS"
sed -n "${START},${END}p" "$MASTER_LIST" > "$SCRATCH/ligand.txt"

if [ ! -s "$SCRATCH/ligand.txt" ]; then
    echo "Task $SLURM_ARRAY_TASK_ID: no ligands in range ${START}-${END} (past end of list) -- nothing to do"
    rmdir "$SCRATCH" 2>/dev/null
    exit 0
fi
echo "Task $SLURM_ARRAY_TASK_ID docking ligands ${START}-${END} ($(wc -l < "$SCRATCH/ligand.txt") ligands)"

#############################################
# Stage inputs into scratch (symlink ligands -- there can be thousands)
cp "$EXECUTABLE" "$CONFIG" "$RECEPTOR" "$SCRATCH/"
while read -r fname; do
    [ -n "$fname" ] && ln -sf "$LIGAND_DIR/$fname" "$SCRATCH/$fname"
done < "$SCRATCH/ligand.txt"

cd "$SCRATCH"

#############################################
# Run
printf "ligand.txt\n" | perl Vina_rigid_P.pl
EXIT_CODE=$?
echo "Finished at: $(date) -- exit code: $EXIT_CODE"

#############################################
# Copy docking results back (skip the symlinked input ligands)
cp *_out.pdbqt "$RESULTS/" 2>/dev/null

#############################################
# Clean up scratch
cd "$SLURM_SUBMIT_DIR"
rm -rf "$SCRATCH"

echo "Results saved to: $RESULTS"
# Exit with $EXIT_CODE directly -- see gsta_docking.sh / docs/implementation_notes.md
exit $EXIT_CODE
