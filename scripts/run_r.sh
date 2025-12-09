#!/bin/bash
#SBATCH --job-name=curve_fitting
#SBATCH --account=jiankang1
#SBATCH --partition=standard
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20g
#SBATCH --output=/home/%u/R/ABMRS/log/output-%x-%A_%a.log
#SBATCH --error=/home/%u/R/ABMRS/log/error-%x-%A_%a.log
#SBATCH --array=1-10

module load R/4.4.0  # Load the R module if needed

SCRIPT_PATH="/home/hejunhui/R/ABMRS/scripts/curve_fitting.r"  # Path to the R script

# Adjust these if you want a different number of replicates per array task
REPS_PER_TASK=${REPS_PER_TASK:-5}

echo "Starting curve_fitting task ${SLURM_ARRAY_TASK_ID} with ${REPS_PER_TASK} replicates"

srun Rscript $SCRIPT_PATH ${SLURM_ARRAY_TASK_ID} ${REPS_PER_TASK}

echo "Done task ${SLURM_ARRAY_TASK_ID} with ${REPS_PER_TASK} replicates"