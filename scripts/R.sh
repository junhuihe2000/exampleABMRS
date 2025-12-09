#!/bin/bash
#SBATCH --job-name=PeakDetection
#SBATCH --account=jiankang1
#SBATCH --mail-user=hejunhui@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-5
#SBATCH --partition=standard
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5g
#SBATCH --output=/home/%u/R/ABMRS/log/output-%x-%j.log
#SBATCH --error=/home/%u/R/ABMRS/log/error-%x-%j.log

module load R/4.4.0  # Load the R module if needed

SCRIPT_PATH="/home/hejunhui/R/ABMRS/fmri/script/peak_detection.R"  # Path to the R script

# Run the R script with the current array index
echo "The ${SLURM_ARRAY_TASK_ID}-th task"
srun Rscript $SCRIPT_PATH ${SLURM_ARRAY_TASK_ID}
echo "The ${SLURM_ARRAY_TASK_ID}-th task is done"