import sys

import model.scenario
from model.scenario import calibrate,foi,art

model.scenario.uid = '2023-02-22'
# model.scenario.N['sam'] = 1000 # DEBUG

# calibrate ---------------------------
# calibrate.run(int(sys.argv[1]))
# calibrate.merge()
# foi ---------------------------------
# calibrate.case = 
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
