#!/bin/bash
#SBATCH --job-name=example
#SBATCH --account=jiankang1
#SBATCH --mail-user=hejunhui@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --partition=standard
#SBATCH --time=10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10g
#SBATCH --output=/home/%u/example/log/output-%x-%j.log
#SBATCH --error=/home/%u/example/log/error-%x-%j.log
#SBATCH --array=1-100

module load python3.11-anaconda/2024.02

echo "The ${SLURM_ARRAY_TASK_ID}-th task"
python example.py --id ${SLURM_ARRAY_TASK_ID}

echo "All done! Finished!"