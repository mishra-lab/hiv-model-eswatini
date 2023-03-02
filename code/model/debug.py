import numpy as np
from model import system,params,target,out,plot,fit,slicers
from utils import dict_list_update
from utils import tarray as ta

t  = system.get_t(tf=2025)

N = 70
T  = target.get_all_esw()
Ps = params.get_n_all(N,seeds=range(N))
Rs = system.run_n(Ps,t,T,para=True)
fit.plot_sets(t,Rs,T)
