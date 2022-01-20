from model import system,params

N  = 10
t  = system.f_t(t1=2025)
Ps = params.get_n_all(N)
Rs = system.run_n(Ps,t,para=False)