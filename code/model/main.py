import sys
from utils import fio,parallel
from model import system,params,target,fit,scenario,out
import numpy as np

# N = 70
# scenario.uid = '2022-01-24'
# t  = system.f_t(t1=2040)
# T1 = target.get_all_esw()
# P1s = params.get_n_all(N,seeds=range(N))
# R1s = system.run_n(P1s,t,T1)
# R1s = system.drop_fails(R1s)[0]
# fit.plot_all(t,R1s,T1,fname=scenario.fname('fig','tmp',N,0)) # DEBUG

scenario.uid = '2022-01-24'
scenario.main(256,sample=False,refit=False)
