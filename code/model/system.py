import numpy as np
from model import params,target,foi
from utils import _,log,deco,ppool,rk4step

def get_t(t0=1980,tf=2025,dt=0.05):
  # define a time vector with some defaults
  return np.round(np.arange(t0,tf+dt,dt),9)

@deco.nowarn
def get_X(X0,t):
  # initialize a population-time array X: all nan but X[0,...] = X0
  X = np.nan * np.ndarray([*t.shape,*X0.shape])
  X[0] = X0
  return X

def run_n(Ps,t=None,T=None,para=True,**kwds):
  # run the model once for each param set P in Ps, usually in parallel
  log(2,'system.run_n: N = '+str(len(Ps)))
  if para:
    fun = lambda P: run(P,t=t,T=T,**kwds)
    return log(-1,ppool().map(fun,Ps))
  else:
    return log(-1,[run(P,t=t,T=T,**kwds) for P in Ps])

def run(P,t=None,T=None,RPts=None,Xk=False):
  # wrapper for solve (below) with some setup & cleanup
  if t is None: t = get_t()
  if RPts is None:
    # RPts are tarray params which we add to R to make some out.functions easier
    RPts = ['PF_condom_t','PF_circum_t','dx_sit','tx_sit','Rtx_ht','unvx_sit','revx_t']
  R = solve(P,t)
  log(3,str(P['id']).rjust(9)+(' ' if R else '!'))
  if not R: # if aborted/failed for any reason
    return {'P':P,'t':t,'ll':-np.inf}
  # sum k (EPA) dimension & possibly remove Xk to save memory
  R['X'] = (R.get('Xk') if Xk else R.pop('Xk')).sum(axis=3)
  R['lls'] = target.get_model_ll(T,R,t,aggr=False) if T else {} # ll per target
  R['ll'] = sum(R['lls'].values()) if T else None # overall ll (log-likelihood)
  if RPts: # evaluate tarray params for all t & reshape so t-dim is first
    R.update({k:np.rollaxis(P[k](t),-1) for k in RPts})
  return R

#@profile
def solve(P,t):
  # solve the model for param set P & time vector t
  X   = get_X(P['X0'],t)               # (t:*, s:2, i:4, k:5, h:6, c:5)
  inc = get_X(np.zeros([4,2,4,2,4]),t) # (t:*, p:4, s:2, i:4, s':2, i':4)
  b_hiv,b_tpaf = True,True # toggles so we only introduce HIV & start TPAF once
  for i in range(1,t.size):
    # Ri = rk4step(X[i-1],t[i-1],(t[i]-t[i-1]),get_dX,P=P)
    Ri = get_dX(X[i-1],t[i-1],P) # DEBUG: Euler
    X[i] = X[i-1] + (t[i] - t[i-1]) * Ri['dX'] # X(t) = X(t-dt) + dt * dX/dt(t-dt)
    inc[i] = Ri['inc']
    if b_hiv and t[i] >= P['t0_hiv']: # introduce HIV
      b_hiv = False
      X[i,:,:,0,:,0] = X[i,:,:,0,0,0,_] * P['PX_h_hiv'][_,_,:]
    if b_tpaf and t[i] >= P['t0_tpaf']: # start accumulating TPAF
      b_tpaf = False
      P['mix_mask'] = P['mix_mask_tpaf']
    if np.any(X[i].sum(axis=2) < 0) or np.any(inc[i] < 0): # abort / fail
      return False
  return {
    'P': P,     # param set
    't': t,     # time vector
    'Xk': X,    # population-time array (keep k dim for now)
    'inc': inc, # incidence-time array
  }

#@profile
def get_dX(X,t,P):
  # initialize dX (fastest)
  dX = 0*X # (s:2, i:4, k:5, h:6, c:5)
  # force of infection - modifies dX internally
  inc = foi.get_apply_inc(dX,X,t,P) # (p:4, s:2, i:4, s':2, i':4)
  # HIV progression
  dXi = X[:,:,:,1:5,0:3] * P['prog_h'] # all hiv & untreated
  dX[:,:,:,1:5,0:3] -= dXi
  dX[:,:,:,2:6,0:3] += dXi
  # CD4 recovery
  dXi = X[:,:,:,3:6,3:5] * P['unprog_h'] # low CD4 & treated
  dX[:,:,:,3:6,3:5] -= dXi
  dX[:,:,:,2:5,3:5] += dXi
  # births & deaths (entry & exit)
  birth, PXe_si, turn = params.solve_turnover(P,t)
  dX[:,:,0,0,0] += X.sum() * birth * PXe_si
  dX -= X * P['death']
  dX -= X * P['death_hc']
  # turnover among activity groups
  dXi = turn[:,:,:,_,_,_] * X[:,:,_,:,:,:]
  dX -= dXi.sum(axis=2) # (s:2, i:4, k:5, h:6, c:5)
  dX += dXi.sum(axis=1) # (s:2, i':4, k:5, h:6, c:5)
  # cascade: diagnosis
  dXi = X[:,:,:,1:6,0] * P['dx_sit'](t) * P['Rdx_scen']
  dX[:,:,:,1:6,0] -= dXi # undiag
  dX[:,:,:,1:6,1] += dXi # diag
  # cascade: treatment
  dXi = X[:,:,:,1:6,1] * P['tx_sit'](t) * P['Rtx_ht'](t) * P['Rtx_scen']
  dX[:,:,:,1:6,1] -= dXi # diag
  dX[:,:,:,1:6,3] += dXi # treat
  # cascade: viral suppression
  dXi = X[:,:,:,1:6,3] * P['vx']
  dX[:,:,:,1:6,3] -= dXi # treat
  dX[:,:,:,1:6,4] += dXi # vls
  # cascade: treatment fail / discontinue
  dXi = X[:,:,:,1:6,4] * P['unvx_sit'](t) * P['Rux_scen']
  dX[:,:,:,1:6,4] -= dXi # vls
  dX[:,:,:,1:6,2] += dXi # fail
  # cascade: viral re-suppression
  dXi = X[:,:,:,1:6,2] * P['revx_t'](t)
  dX[:,:,:,1:6,2] -= dXi # fail
  dX[:,:,:,1:6,4] += dXi # vls
  return {
    'dX': dX,
    'inc': inc,
  }
