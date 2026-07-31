#!/bin/bash -l
#SBATCH --account=pawsey1376
#SBATCH --job-name=gst_py_16
#SBATCH --partition=work
#SBATCH --nodes=1
##SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00

# Initialise Conda
source $MYSOFTWARE/miniconda3/etc/profile.d/conda.sh

# Activate your environment
conda activate gst_ml

# --- OpenMP settings ---
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# --- MPI settings ---
export MPICH_OFI_STARTUP_CONNECT=1
# export MPICH_OFI_VERBOSE=1        # Uncomment for debugging
export MPICH_ENV_DISPLAY=1        # Useful for testing hybrid setup
# export MPICH_MEMORY_REPORT=1

# --- Run ---
echo "Running on $SLURM_NNODES nodes with $SLURM_NTASKS MPI ranks × $OMP_NUM_THREADS threads each."

echo "Node(s):         $SLURM_JOB_NUM_NODES"
echo "Task(s):         $SLURM_NTASKS"
echo "Threads/task:    $OMP_NUM_THREADS"
echo "Python:          $(which python)"
echo "Environment:     $CONDA_DEFAULT_ENV"
 
srun -N $SLURM_NNODES -n $SLURM_NTASKS -c $OMP_NUM_THREADS python pipeline_HAL20260712_ensemble.py
