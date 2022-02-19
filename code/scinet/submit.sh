#!/bin/bash
#SBATCH --job-name hiv-esw
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=1:00:00

module load python/3.8
# pip install pathos PyPDF2
export PYTHONPATH=.:$PYTHONPATH
python3.8 scinet/main.py $SLURM_ARRAY_TASK_ID

# using this script on scinet:
# sbatch -a 0-10 submit.sh
