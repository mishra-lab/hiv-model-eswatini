import sys,os
os.environ['MPLCONFIGDIR'] = sys.argv[1]

from utils import parallel
from model import scenario

parallel.cpus = 80
scenario.uid = '2022-02-16'
scenario.N['size']  = 10000
scenario.N['batch'] = int(sys.argv[1])
scenario.main_fit(cf=False)