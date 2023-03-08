import sys,os
os.environ['MPLCONFIGDIR'] = os.path.join('.tmp',sys.argv[1])

from utils import parallel
import model.scenario
from model.scenario import calibrate,tpaf,foi,art

parallel.cpus = 80
model.scenario.uid = '2023-03-03'
# model.scenario.N['sam'] = 1000 # DEBUG

# calibrate ---------------------------
# calibrate.run(int(sys.argv[1]))
# calibrate.merge()
# calibrate.rerun()
# tpaf --------------------------------
# tpaf.run(sys.argv[2])
# foi ---------------------------------
# calibrate.case = sys.argv[2]
# calibrate.run(int(sys.argv[1]),foi_mode=calibrate.case)
# calibrate.merge()
# foi.run_ep()
# foi.run_fit(case)
# foi.run_tpaf(case)
# art ---------------------------------
# for case in art.cases:
#   art.run_refit(case,int(sys.argv[1]))
# art.merge_refit()
# art.rerun_refit()
# art.run_sens()
