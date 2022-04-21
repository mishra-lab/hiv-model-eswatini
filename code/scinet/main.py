import sys,os
os.environ['MPLCONFIGDIR'] = os.path.join('.tmp',sys.argv[1])

from utils import parallel
import model.scenario
from model.scenario import calibrate

parallel.cpus = 80
model.scenario.uid = '2022-04-20'
model.scenario.N['b'] = int(sys.argv[1])

calibrate.run()
# calibrate.merge()
# calibrate.rerun()