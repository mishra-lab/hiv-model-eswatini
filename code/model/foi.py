import numpy as np
from utils import _,deco,linear_comb

tol = 1e-7
modes = [
  'lin', # no consideration of partnership duration
  'bpd', # binomial per-partnership-duration
  'bpy', # binomial per-partnership-year
  'bmd', # binomial concurrent partnership-duration
  'bmy', # binomial concurrent partnership-year
  'fpe', # force of infection partnership exclusion
]
mode = 'fpe'

#@profile
def get_beta(P,t):
  # return.shape = (a:2, p:4, s:2, i:4, s':2, i':4, h':6, c':5)
  RbA_condom = linear_comb(P['PF_condom_t'](t) * P['RPF_condom_a'], P['Rbeta_condom'], 1)
  RbA_circum = linear_comb(P['PF_circum_t'](t), P['Rbeta_circum'], 1)
  P_gud_t = P['P_gud_t'](t) * P['P_gud']
  Rbeta_gud_sus = linear_comb(P_gud_t,P['Rbeta_gud_sus'],1).reshape([1,1,2,4,1,1,1,1])
  Rbeta_gud_inf = linear_comb(P_gud_t,P['Rbeta_gud_inf'],1).reshape([1,1,1,1,2,4,1,1])
  return P['beta_a'] * Rbeta_gud_sus * Rbeta_gud_inf * RbA_condom * RbA_circum

@deco.nowarn
#@profile
def get_mix(XC,P):
  # return.shape = (p:4, s:2, i:4, s':2, i':4)
  M0 = XC[:,0,:,_] * XC[:,1,_,:] / XC.sum(axis=2).mean(axis=1)[:,_,_] + tol/10
  m1 = M0.sum(axis=1)
  m2 = M0.sum(axis=2)
  # print(m1 / XC[:,1,:]) # DEBUG == 1, unless XC unbalanced
  # print(m2 / XC[:,0,:]) # DEBUG == 1, unless XC unbalanced
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
  # n.b. returning on population-level scale
  return P['mix']

@deco.nowarn
#@profile
def get_apply_inc(dX,X,t,P):
  # return.shape = (p:4, s:2, i:4, s':2, i':4)
  # NOTE: if mode in ['lin','fpe','bpd','bpy']: return absolute infections
  #       if mode in ['bmd','bmy']: return *probability* of infection (aggr must be deferred)
  # define partner numbers (C) / rates (Q) for mixing, + total acts (A) for binomial models
  if mode in ['fpe']: # C
    C_psik = P['C_psi'] - P['aC_pk']
  elif mode in ['lin','bpy','bmy']: # C, Q, Q
    C_psik = P['C_psi']
    A_ap = P['F_ap']
  elif mode in ['bpd','bmd']: # Q, Q
    C_psik = P['C_psi'] / P['dur_p'][:,_,_,_]
    A_ap = P['F_ap'] * P['dur_p'][_,:]
  # compute the mixing
  XC = (X[_,:,:,:,:,:] * C_psik[:,:,:,:,_,_]).sum(axis=3)
  SXC = XC.sum(axis=(3,4))  # shape = (p:4, s:2, i:4)
  PXC_hc = XC / (SXC[:,:,:,_,_] + tol/10) # shape = (p:4, s:2, i:4, h:6, c:5)
  mix = get_mix(SXC,P) * PXC_hc[:,:,:,0,0,_,_] * P['mix_mask']
  # compute per-act probability
  beta = get_beta(P,t)
  # force of infection: linear vs binomial
  if mode in ['lin','fpe']:
    Fbeta = np.sum(beta * P['F_ap'][:,:,_,_,_,_,_,_], axis=0)
    inc = Fbeta * mix[:,:,:,:,:,_,_] * PXC_hc[:,_,_,:,:,:,:] # absolute infections
  elif mode in ['bpd','bpy','bmd','bmy']:
    Abeta = 1 - np.prod((1 - beta) ** A_ap[:,:,_,_,_,_,_,_], axis=0)
    XAbeta = np.sum(Abeta * PXC_hc[:,_,_,:,:,:,:], axis=(5,6))
    if mode in ['bpd','bpy']:
      inc = XAbeta * mix # absolute infections
    elif mode in ['bmd','bmy']:
      inc = 1 - (1 - XAbeta) ** (mix / X[_,:,:,0,0,0,_,_]) # probability of infection
  # aggregating & applying to dX
  if mode in ['fpe']:
    dXi = inc.sum(axis=(3,4,5,6)) # acquisition: (p:4, s:2, i:4)
    dX[:,:,0 ,0,0] -= dXi.sum(axis=0)
    dX[:,:,1:,1,0] += np.moveaxis(dXi,0,2)
    dXi = inc.sum(axis=(1,2)) # transmission: (p:4, s':2, i':4, h':6, c':5)
    dX[:,:,0 ,:,:] -= dXi.sum(axis=0)
    dX[:,:,1:,:,:] += np.moveaxis(dXi,0,2)
    dXi = X[:,:,1:,:,:] / P['dur_p'][_,_,:,_,_] # new partnerships: (s:2, i:4, k:4, h:6, c:5)
    dX[:,:,1:,:,:] -= dXi
    dX[:,:,0 ,:,:] += dXi.sum(axis=2)
    return inc.sum(axis=(5,6))
  # all other models don't use dimension "k" in X
  elif mode in ['lin']:
    inc = inc.sum(axis=(5,6))
    dXi = aggr_inc(inc,axis=(0,3,4))
  elif mode in ['bpd','bpy']:
    dXi = aggr_inc(inc,axis=(0,3,4))
  elif mode in ['bmd','bmy']:
    dXi = aggr_inc(inc,axis=(0,3,4),Xsus=X[:,:,0,0,0])
  dX[:,:,0,0,0] -= dXi
  dX[:,:,0,1,0] += dXi
  return inc

#@profile
def aggr_inc(inc,axis,Xsus=1.,keepdims=False):
  # returns absolute infections after appropriately aggregating "inc"
  # Xsus only needed if mode in ['bmy','bmd']
  if mode in ['lin','fpe','bpy','bpd']:
    return inc.sum(axis=axis,keepdims=keepdims)
  if mode in ['bmy','bmd']:
    return (1 - np.prod(1 - inc,axis=axis,keepdims=keepdims)) * Xsus
