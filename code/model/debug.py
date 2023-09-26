import numpy as np
from model import system,params,target,out,plot,fit,slicers
from utils import dict_list_update
from utils import tarray as ta

t  = system.get_t(tf=2025)
T  = target.get_all_esw()
Ps = params.get_n_all(range(70))
Rs = system.run_n(Ps,t,T,para=True)
fit.plot_sets(t,Rs,T=T,debug=True)
