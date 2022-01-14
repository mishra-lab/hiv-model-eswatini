import sys,os
UID = sys.argv[1]
N   = int(sys.argv[2])
B   = int(sys.argv[3])
os.environ['MPLCONFIGDIR'] = str(B)

from utils import parallel
from model import scenario,handfit

parallel.cpus = 80
handfit.tfname = str(B)+'/'+handfit.tfname
scenario.main(N=N,N0=N*B,cf=False)
