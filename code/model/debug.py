import numpy as np
from model import target,system,params,fit

N = 10
T = target.get_all_esw()
t  = system.get_t(tf=2021)
Ps = params.get_n_all(N,seeds=range(N))
Rs = system.run_n(Ps,t,T,para=True)
fit.plot_debug(t,Rs,T)