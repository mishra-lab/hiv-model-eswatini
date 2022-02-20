import sys,os
os.environ['MPLCONFIGDIR'] = os.path.join('.tmp',sys.argv[1])

from utils import parallel
from model import scenario

parallel.cpus = 80
scenario.uid = '2022-02-20'
scenario.N['b'] = int(sys.argv[1])

# scenario.run_calibrate()
# scenario.merge_calibrate()
scenario.run_refit()
# scenario.merge_refit()
# scenario.expo_refit()
# scenario.run_sens()
