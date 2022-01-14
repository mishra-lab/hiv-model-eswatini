#!/bin/bash
#SBATCH --job-name hiv-eswatini
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=0:30:00

N=2500
module load python/3.8
pip install pathos
export PYTHONPATH=.:$PYTHONPATH
python3.8 scinet/main.py $N $SLURM_ARRAY_TASK_ID

# using this script on scinet:
# sbatch -a 0-10 submit.sh