import numpy as np
from model import foi,target
from utils import _,deco,parallel,log,rk4step,nan_to_value

def get_t(t0=1980,tf=2050,dt=0.05):
  return np.round(np.arange(t0,tf+dt,dt),9)

@deco.nowarn
def get_X(X0,t):
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
  if t is None: t = get_t()
  if RPts is None: RPts = ['PF_condom_t','PF_circum_t','P_gud_t','dx_t','tx_t','Rtx_ht']
  R = solve(P,t)
  log(3,str(P['seed']).rjust(6)+(' ' if R else '!'))
  if not R:
    return R
  if T is not None:
    R['ll'] = target.get_model_ll(T,R,t,interval=interval)
    R['P']['ll'] = R['ll']
  if RPts:
    for RPt in RPts:
      R[RPt] = np.rollaxis(P[RPt](t),-1)
  R['P'].pop('inc')
  R['P'].pop('beta')
  return R

def solve(P,t):
  X   = get_X(P['X0'],t)
  inc = get_X(np.zeros([4,2,4,2,4]),t)
  t0_hiv = int(P['t0_hiv']) # TODO: avoid int here
  t0_tpaf = int(P['t0_tpaf'])
  for i in range(1,t.size):
    # Ri = rk4step(X[i-1],t[i-1],(t[i]-t[i-1]),get_dX,P=P)
    Ri = get_dX(X[i-1],t[i-1],P) # DEBUG: Euler
    X[i] = X[i-1] + (t[i] - t[i-1]) * Ri['dX']
    inc[i] = Ri['inc']
    if t[i] == t0_hiv: # introduce HIV
      X[i,:,:,0,:,0] = X[i,:,:,0,0,0,_] * P['PX_h_hiv'][_,_,:]
    if t[i] == t0_tpaf: # start accumulating tPAF
      P['mix_mask'] = P['mix_mask_tpaf']
    if np.any(X[i].sum(axis=2) < 0): # abort / fail
      return False
    # check.all(P,X[i],t[i])
  return {
    'P': P,
    'X': X.sum(axis=3), # sum_k
    't': t,
    'inc': inc,
  }

#@profile
def get_dX(X,t,P):
  # initialize
  dX = 0*X
  # force of infection
  # ([a:2], p:4, s:2, i:4, s':2, i':4, h':6, c':5)
  P['beta'] = foi.get_beta(P,t)
  P['inc']  = foi.get_inc(P,X,t) # TODO: * P['mix_mask']
  foi.apply_inc(dX,P['inc'],X)
  # forming new partnerships [foi.mode = 'fpe' only]
  dXi = X[:,:,1:,:,:] / P['dur_p'][_,_,:,_,_]
  dX[:,:,1:,:,:] -= dXi
  dX[:,:,0 ,:,:] += dXi.sum(axis=2)
  # HIV transitions
  dXi = X[:,:,:,1:5,0:3] * P['prog_h'] # all hiv & untreated
  dX[:,:,:,1:5,0:3] -= dXi
  dX[:,:,:,2:6,0:3] += dXi
  # CD4 recovery
  dXi = X[:,:,:,3:6,3:5] * P['unprog_h'] # low CD4 & treated
  dX[:,:,:,3:6,3:5] -= dXi
  dX[:,:,:,2:5,3:5] += dXi
  # births & deaths
  dX[:,:,0,0,0] += X.sum() * P['PX_si'] * P['birth_t'](t)
  dX -= X * P['death']
  dX -= X * P['death_hc']
  # turnover
  dXi = get_turnover(P,X)
  dX -= dXi.sum(axis=2) # (s:2, i:4, k:5, h:6, c:5)
  dX += dXi.sum(axis=1) # (s:2, i':4, k:5, h:6, c:5)
  # cascade: diagnosis
  dXi = X[:,:,:,1:6,0] * P['dx_t'](t) * P['Rdx_si'] * P['Rdx_scen']
  dX[:,:,:,1:6,0] -= dXi # undiag
  dX[:,:,:,1:6,1] += dXi # diag
  # cascade: treatment
  dXi = X[:,:,:,1:6,1] * P['tx_t'](t) * P['Rtx_ht'](t) * P['Rtx_si'] * P['Rtx_scen']
  dX[:,:,:,1:6,1] -= dXi # diag
  dX[:,:,:,1:6,3] += dXi # treat
  # cascade: VLS
  dXi = X[:,:,:,1:6,3] * P['vx']
  dX[:,:,:,1:6,3] -= dXi # treat
  dX[:,:,:,1:6,4] += dXi # vls
  # cascade: unlink
  dXi = X[:,:,:,1:6,4] * P['unvx_t'](t) * P['Rux_scen']
  dX[:,:,:,1:6,4] -= dXi # vls
  dX[:,:,:,1:6,2] += dXi # unlink
  # cascade: relink
  dXi = X[:,:,:,1:6,2] * P['retx']
  dX[:,:,:,1:6,2] -= dXi # unlink
  dX[:,:,:,1:6,4] += dXi # vls
  return {
    'dX': dX,
    'inc': foi.aggr_inc(P['inc'],axis=(5,6)),
  }

@deco.nowarn
#@profile
def get_turnover(P,X):
  # return.shape = (s:2, i:4, i':4, k:5, h:6, c:5)
  turn = P['turn_sii'][:,:,:,_,_,_] * X[:,:,_,:,:,:]
  # P['ORturn_sus:hiv'] = 0 # DEBUG
  if np.any(X[:,:,:,1:,:]): # HIV introduced
    Xhiv = X[:,:,:,1:,:].sum(axis=(2,3,4))
    # odds of turnover aomng sus vs hiv (source-group-specific)
    Osus = P['ORturn_sus:hiv'] * nan_to_value(X[:,:,0,0,0] / Xhiv, 1)[:,:,_]
    Phc_hiv = nan_to_value(X[:,:,:,1:,:] / Xhiv[:,:,_,_,_], 0)
    turn_hiv = turn.sum(axis=(3,4,5)) / (1 + Osus)
    turn[:,:,:,:,1:,:] = turn_hiv[:,:,:,_,_,_] * Phc_hiv[:,:,_,:,:,:]
    turn[:,:,:,0,0,0]  = turn_hiv * Osus
  return turn
