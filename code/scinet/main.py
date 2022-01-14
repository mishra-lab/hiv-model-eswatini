import sys
from utils import parallel
from model import scenario

N  = sys.argv[1]
N0 = sys.argv[2] * N
parallel.cpus = 40
scenario.uid = '2022-01-14'
scenario.main(N=N,N0=N0,cf=False)
