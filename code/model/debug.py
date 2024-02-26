import numpy as np
from model import system,params,target,fit

t  = system.get_t(tf=2025)
T  = target.get_all_esw()
Ps = params.get_n_all(range(70))
Rs = system.run_n(Ps,t,T,para=True,Xk=True)
fit.plot_sets(t,Rs,T=T,debug=True)
