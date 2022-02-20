import sys
from utils import parallel
from model import scenario
import numpy as np

scenario.uid = '2022-02-20'
# scenario.N['cal'] = 1000 # DEBUG
scenario.N['b'] = int(sys.argv[1])

scenario.run_calibrate()
scenario.merge_calibrate()
scenario.run_refit()
scenario.merge_refit()
scenario.expo_refit()
scenario.run_sens()