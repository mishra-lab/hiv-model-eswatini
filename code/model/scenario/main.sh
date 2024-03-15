#!/bin/bash
#SBATCH --job-name hiv-esw
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=0:15:00
#SBATCH --array=0

# setup once {scinet} ----------------------------------------------------------
# pip install --upgrade scipy pathos PyPDF2

# setup every time -------------------------------------------------------------
module load python/3.8 # {scinet}
PYCODE="""
from model.scenario import imis,art,foi
b = $SLURM_ARRAY_TASK_ID # {scinet}
"""
export PYTHONPATH=.:$PYTHONPATH
export MPLCONFIGDIR=.tmp # {scinet}

# base calibration -------------------------------------------------------------
# PYCODE+="imis.run(case='base',b=b)"                 # 0:45:00; array=0-99
# PYCODE+="imis.sample_post(case='base')"             # 0:15:00
# PYCODE+="imis.rerun(case='base')"                   # 0:15:00

# art calibration & scenarios --------------------------------------------------
# PYCODE+="art.run_rf(b=b)"                           # 1:30:00; array=0-3
# PYCODE+="art.rerun_rf()"                            # 0:45:00
# PYCODE+="art.run_ss()"                              # 0:45:00

# foi calibration & scenarios --------------------------------------------------
# PYCODE+="case = 'foi-TODO'; mode = case.replace('foi-','');"
# PYCODE+="foi.run_ep()"                              # 0:30:00
# PYCODE+="imis.run(case=case,b=b,foi_mode=mode)"     # 0:45:00; array=0-99
# PYCODE+="imis.sample_post(case=case)"               # 0:15:00
# PYCODE+="imis.rerun(case=case)"                     # 0:15:00
# PYCODE+="foi.run_tpaf(case=case)"                   # 4:00:00

# run every time ---------------------------------------------------------------
echo "$PYCODE"
python3 -c "$PYCODE"
