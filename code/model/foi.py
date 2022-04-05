import numpy as np
from utils import _,deco,linear_comb,nan_to_value

# TODO: merge foi.py with system.py

#@profile
def f_lambda_pp(P,t):
  # return.shape = (p:4, s:2, i:4, s':2, i':4, h':6, c':5)
  # compute some stratified effect modifiers
  RbA_condom = linear_comb(P['PF_condom_t'](t) * P['RPF_condom_a'], P['Rbeta_condom'], 1)
  RbA_circum = linear_comb(P['PF_circum_t'](t), P['Rbeta_circum'], 1)
  P_gud_t   = P['P_gud_t'](t) * P['P_gud']
  Rbeta_gud_sus = linear_comb(P_gud_t,P['Rbeta_gud_sus'],1).reshape([1,1,2,4,1,1,1,1])
  Rbeta_gud_inf = linear_comb(P_gud_t,P['Rbeta_gud_inf'],1).reshape([1,1,1,1,2,4,1,1])
  beta_a = P['beta_a'] * Rbeta_gud_sus * Rbeta_gud_inf
  # per-act transm prob <= .5
  beta_a = np.minimum(beta_a,.5)
  # multiply all stratified elements together
  return np.sum(beta_a * P['F_ap'] * RbA_condom * RbA_circum, axis=0)

@deco.nowarn
#@profile
def f_lambda(P,X):
  # return.shape = (p:4, s:2, i:4, s':2, i':4, h':6, c':5)
  tol = 1e-7
  # total partners offered
  XC = (X[_,:,:,:,:,:] * P['C_psik'][:,:,:,:,_,_]).sum(axis=3)
  XCs = XC.sum(axis=(3,4))
  PXC_hc = XC / (XCs[:,:,:,_,_] + tol/10)
  # random mixing: shape = (p:4, s:2, i:4, i':4)
  M0 = XCs[:,0,:,_] * XCs[:,1,_,:] / XCs.sum(axis=2).mean(axis=1)[:,_,_] + tol/10
  m1 = M0.sum(axis=1)
  m2 = M0.sum(axis=2)
  # print(m1 / XC[:,1,:]) # DEBUG == 1, unless XC unbalanced
  # print(m2 / XC[:,0,:]) # DEBUG == 1, unless XC unbalanced
  # preferences + margin fix (iterative proportional fitting)
  M = M0 * np.exp(P['pref_pii'])
  for k in range(100):
    r1 = m1 / M.sum(axis=1)
    M *= r1[:,_,:]
    r2 = m2 / M.sum(axis=2)
    M *= r2[:,:,_]
    if (abs(r1-1) < tol).all() and (abs(r2-1) < tol).all():
      break
  M[abs(M)<tol] = 0
  P['mix'][:,0,:,1,:] = M
  P['mix'][:,1,:,0,:] = M.swapaxes(1,2)
  # per-partner rate * mixing * susceptible * infectious
  return P['lambda_pp'] * P['mix'][:,:,:,:,:,_,_] * PXC_hc[:,:,:,0,0,_,_,_,_] * PXC_hc[:,_,_,:,:,:,:]

@deco.nowarn
#@profile
def f_turnover(P,X):
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
