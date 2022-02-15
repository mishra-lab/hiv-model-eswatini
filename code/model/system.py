import numpy as np
from model import foi,target
from utils import _,deco,parallel

def f_t(t0=1980,t1=2050,dt=0.1):
  return np.round(np.arange(t0,t1+dt,dt),9)

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
  if para:
    fun = lambda P: run(P,t=t,T=T,**kwds)
    return parallel.ppool(len(Ps)).map(fun,Ps)
  else:
    return [run(P,t=t,T=T,**kwds) for P in Ps]
  print(flush=True)

def run(P,t=None,T=None,RPts=None,interval=None):
  if t is None: t = f_t()
  if RPts is None: RPts = ['PA_condom_t','PA_circum_t','P_gud_t','dx_t','tx_t','Rtx_ht']
  print(str(P['seed']).rjust(6),end=' ',flush=True)
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
  esc = f_X(np.ones([4,2,4,2,4]),t)
  t0_hiv = int(P['t0_hiv'])
  for i in range(1,t.size):
    # TODO: use array.dtfun?
    R = f_dX(P,X[i-1],t[i-1])
    X[i] = X[i-1] + (t[i] - t[i-1]) * R['dX']
    esc[i] = R['esc']
    if t[i] == t0_hiv:
      X[i] = X[i,:,:,_,0,:] * P['PX_h_hiv']
    if np.any(X[i] < 0): # abort / fail
      return False
    # check.all(P,X[i],t[i])
  return {
    'P': P,
    'X': X,
    't': t,
    'esc': esc,
  }

#@profile
def f_dX(P,X,t):
  P['beta_p']  = foi.f_beta_p(P,t)
  P['beta_pp'] = foi.f_beta_pp(P,X)
  P['mix']     = foi.f_mix(P,X)
  # initialize
  dX = 0*X
  # force of infection
  esc = foi.f_lambda(P,X)
  inc = 1 - np.prod(esc,axis=(0,3,4))
  dXi = X[:,:,0,0] * inc
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
  dX[:,:,_,0,_,0] += X.sum() * P['birth_si']
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
  dX[:,:,1:6,3] += dXi # treat
  return {
    'dX': dX,
    'esc': esc,
  }
