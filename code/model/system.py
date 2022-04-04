import numpy as np
from model import foi,target
from utils import _,rk4step,deco,parallel,log

def f_t(t0=1980,tf=2050,dt=0.1):
  return np.round(np.arange(t0,tf+dt,dt),9)

@deco.nowarn
def f_X(X0,t):
  X = np.nan * np.ndarray([*t.shape,*X0.shape])
  X[0] = X0
  return X

def drop_fails(*Rss):
  # Rss is a list of lists of R, which are assumed to be paired.
  # e.g. Rss[0][i] is fit i from scenario A, Rss[1][i] is (paired) fit i from scenario B
  # We drop all fits that have any fail across scenarios
  oks = [all(Ris) for Ris in zip(*Rss)]
  return tuple([R for (R,ok) in zip(Rs,oks) if ok] for Rs in Rss)

def run_n(Ps,t=None,T=None,para=True,**kwds):
  log(2,'system.run_n: '+str(len(Ps)))
  if para:
    fun = lambda P: run(P,t=t,T=T,**kwds)
    Rs = parallel.ppool(len(Ps)).map(fun,Ps)
  else:
    Rs = [run(P,t=t,T=T,**kwds) for P in Ps]
  log(1)
  return Rs

def run(P,t=None,T=None,RPts=None,interval=None):
  if t is None: t = f_t()
  if RPts is None: RPts = ['PF_condom_t','PF_circum_t','P_gud_t','dx_t','tx_t','Rtx_ht']
  log(3,str(P['seed']).rjust(6))
  R = solve(P,t)
  if not R:
    return R
  if T is not None:
    R['ll'] = target.get_model_ll(T,R,t,interval=interval)
    R['P']['ll'] = R['ll']
  if RPts:
    for RPt in RPts:
      R[RPt] = np.rollaxis(P[RPt](t),-1)
  return R

def solve(P,t):
  X   = f_X(P['X0'],t)
  inc = f_X(np.zeros([4,2,4,2,4]),t)
  t0_hiv = int(P['t0_hiv'])
  t0_tpaf = int(P['t0_tpaf'])
  for i in range(1,t.size):
    Ri = rk4step(X[i-1],t[i-1],(t[i]-t[i-1]),f_dX,P=P)
    # Ri = f_dX(X[i-1],t[i-1],P) # DEBUG: Euler
    X[i] = X[i-1] + (t[i] - t[i-1]) * Ri['dX']
    inc[i] = Ri['inc']
    if t[i] == t0_hiv: # introduce HIV
      X[i] = X[i,:,:,_,0,:] * P['PX_h_hiv']
    if t[i] == t0_tpaf: # start accumulating tPAF
      P['mix_mask'] = P['mix_mask_tpaf']
    if np.any(X[i] < 0): # abort / fail
      return False
    # check.all(P,X[i],t[i])
  return {
    'P': P,
    'X': X,
    't': t,
    'inc': inc,
  }

#@profile
def f_dX(X,t,P):
  P['lambda_p'] = foi.f_lambda_p(P,t)
  P['mix']  = foi.f_mix(P,X) * P['mix_mask']
  # initialize
  dX = 0*X
  # force of infection
  inc = foi.f_lambda(P,X)
  dXi = X[:,:,0,0] * inc.sum(axis=(0,3,4))
  dX[:,:,0,0] -= dXi # sus
  dX[:,:,1,0] += dXi # acute undiag
  # HIV transitions
  dXi = X[:,:,1:5,0:3] * P['prog_h'] # all hiv & untreated
  dX[:,:,1:5,0:3] -= dXi
  dX[:,:,2:6,0:3] += dXi
  # CD4 recovery
  dXi = X[:,:,3:6,3:5] * P['unprog_h']
  dX[:,:,3:6,3:5] -= dXi
  dX[:,:,2:5,3:5] += dXi
  # births & deaths
  dX[:,:,0,0] += X.sum() * P['PX_si'] * P['birth_t'](t)
  dX -= X * P['death']
  dX -= X * P['death_hc']
  # turnover
  dXi = foi.f_turnover(P,X)
  dX -= dXi.sum(axis=2)
  dX += dXi.sum(axis=1)
  # cascade: diagnosis
  dXi = X[:,:,1:6,0] * P['dx_t'](t) * P['Rdx_si'] * P['Rdx_scen']
  dX[:,:,1:6,0] -= dXi # undiag
  dX[:,:,1:6,1] += dXi # diag
  # cascade: treatment
  dXi = X[:,:,1:6,1] * P['tx_t'](t) * P['Rtx_ht'](t) * P['Rtx_scen']
  dX[:,:,1:6,1] -= dXi # diag
  dX[:,:,1:6,3] += dXi # treat
  # cascade: VLS
  dXi = X[:,:,1:6,3] * P['vx']
  dX[:,:,1:6,3] -= dXi # treat
  dX[:,:,1:6,4] += dXi # vls
  # cascade: unlink
  dXi = X[:,:,1:6,4] * P['unvx_t'](t) * P['Rux_scen']
  dX[:,:,1:6,4] -= dXi # vls
  dX[:,:,1:6,2] += dXi # unlink
  # cascade: relink
  dXi = X[:,:,1:6,2] * P['retx']
  dX[:,:,1:6,2] -= dXi # unlink
  dX[:,:,1:6,4] += dXi # vls
  return {
    'dX': dX,
    'inc': inc,
  }
