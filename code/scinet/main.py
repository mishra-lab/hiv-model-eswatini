import sys,os
os.environ['MPLCONFIGDIR'] = os.path.join('.tmp',sys.argv[1])

from utils import parallel
import model.scenario
from model.scenario import calibrate,foi,art

parallel.cpus = 80
model.scenario.uid = '2022-04-20'
model.scenario.N['b'] = int(sys.argv[1])

# calibrate ---------------------------
# calibrate.run()
# calibrate.merge()
# calibrate.rerun()
# foi ---------------------------------
# case = ['bpd','bpy','bmy','base'][int(sys.argv[1])]
# foi.run_ep()
# foi.run_fit(case)
# foi.run_tpaf(case)
# art ---------------------------------
case = art.cases[0]
art.run_refit(case)
# art.rerun_refit()