from model import system,params

N = 10
t  = system.f_t(tf=2021)
Ps = params.get_n_all(N)
Rs = system.run_n(Ps,t,para=False)