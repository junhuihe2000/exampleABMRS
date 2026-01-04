#!/bin/bash
#SBATCH --job-name=sfab
#SBATCH --account=jiankang1
#SBATCH --partition=standard
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20g
#SBATCH --output=/home/%u/R/ABMRS/log/output-%x-%A_%a.log
#SBATCH --error=/home/%u/R/ABMRS/log/error-%x-%A_%a.log
#SBATCH --array=1-10

module load R/4.4.0  # Load the R module if needed

# include site library and your home library
export R_LIBS_BASE="/sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/library"   # use the path you saw in step 1
export R_LIBS_SITE="/sw/pkgs/arc/stacks/gcc/13.2.0/Rtidyverse/4.4.0/2024-05-11"   # use the path you saw in step 1
export R_LIBS_USER="/home/hejunhui/R/x86_64-pc-linux-gnu-library/4.4"                           # your personal installs
export R_LIBS="$R_LIBS_USER:$R_LIBS_SITE:$R_LIBS_BASE"

SCRIPT_PATH="/home/hejunhui/R/ABMRS/scripts/surface_fitting_appendix_bmars.r"  # Path to the R script

# Adjust these if you want a different number of replicates per array task
REPS_PER_TASK=${REPS_PER_TASK:-5}

echo "Starting surface_fitting_appendix_bmars task ${SLURM_ARRAY_TASK_ID} with ${REPS_PER_TASK} replicates"

srun Rscript $SCRIPT_PATH ${SLURM_ARRAY_TASK_ID} ${REPS_PER_TASK}

echo "Done task ${SLURM_ARRAY_TASK_ID} with ${REPS_PER_TASK} replicates"