from model import system,params

N = 10
t  = system.get_t(tf=2025)
Ps = params.get_n_all(range(N))
Rs = system.run_n(Ps,t,para=False)
