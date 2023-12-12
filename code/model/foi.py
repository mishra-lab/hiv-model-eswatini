import numpy as np
from utils import _,deco,linear_comb
from model import tol

foi_modes = [
  'lin', # linear
  'rd', # rate duration
  'ry', # rate year
  'py', # proportion year
  'base', # proposed model (epa)
]

#@profile
def get_beta(P,t):
  # return.shape = (a:2, p:4, s:2, i:4, s':2, i':4, h':6, c':5)
  RbA_condom = linear_comb(P['PF_condom_t'](t) * P['RPF_condom_a'], P['Rbeta_condom'], 1)
  RbA_circum = linear_comb(P['PF_circum_t'](t), P['Rbeta_circum'], 1)
  P_gud_t = P['P_gud'] * P['RP_gud_t'](t)
  Rbeta_gud_sus = linear_comb(P_gud_t,1+P['aRbeta_gud_sus'],1).reshape([1,1,2,4,1,1,1,1])
  Rbeta_gud_inf = linear_comb(P_gud_t,1+P['aRbeta_gud_inf'],1).reshape([1,1,1,1,2,4,1,1])
  return np.minimum(.5, P['beta_a'] * Rbeta_gud_sus * Rbeta_gud_inf * RbA_condom * RbA_circum)

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
  # NOTE: if foi_mode in ['base','li','rd','ry']: return *absolute* infections (not per susceptible)
  #       if foi_mode in ['py']: return *probability* of infection (aggr must be deferred)
  # partner numbers (K) or rates (Q)
  if P['foi_mode'] in ['base']: # K
    C_psik = P['K_psi'] - P['aK_pk']
    fix_XK(X,P)
  elif P['foi_mode'] in ['lin']: # K
    C_psik = P['K_psi']
  elif P['foi_mode'] in ['rd']: # Q
    C_psik = P['K_psi'] / P['dur_p'][:,_,_,_]
    A_ap = P['F_ap'] * P['dur_p']
  elif P['foi_mode'] in ['ry','py']: # Q_1
    C_psik = P['K_psi'] / P['dur_p_1'][:,_,_,_]
    A_ap = P['F_ap'] * P['dur_p_1']
  # compute mixing
  XC = (X[_,:,:,:,:,:] * C_psik[:,:,:,:,_,_]).sum(axis=3)
  XC[XC<0] = 0
  SXC = XC.sum(axis=(3,4))  # shape = (p:4, s:2, i:4)
  PXC_hc = XC / (SXC[:,:,:,_,_] + tol/10) # shape = (p:4, s:2, i:4, h:6, c:5)
  # note: PXC_hc = PX_hc, i.e. X / X.sum(axis=(2,3)), unless foi_mode = 'base'
  mix = get_mix(SXC,P) * P['mix_mask']
  # compute per-act probability
  beta = get_beta(P,t) # shape: (a:2, p:4, s:2, i:4, s':2, i':4, h':6, c':5)
  # compute & apply force of infection
  if P['foi_mode'] in ['base']:
    Fbeta = (beta * P['F_ap'][:,:,_,_,_,_,_,_]).sum(axis=0)
    inc = mix[:,:,:,:,:,_,_] * PXC_hc[:,:,:,0,0,_,_,_,_] * Fbeta * PXC_hc[:,_,_,:,:,:,:]
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
  elif P['foi_mode'] in ['lin']:
    Fbeta = (beta * P['F_ap'][:,:,_,_,_,_,_,_]).sum(axis=0)
    inc = mix * PXC_hc[:,:,:,0,0,_,_] * (Fbeta * PXC_hc[:,_,_,:,:,:,:]).sum(axis=(5,6))
    dXi = aggr_inc(inc,P['foi_mode'],axis=(0,3,4))
  elif P['foi_mode'] in ['rd','ry']:
    B_p = 1 - ((1 - beta) ** A_ap[:,:,_,_,_,_,_,_]).prod(axis=0)
    inc = mix * PXC_hc[:,:,:,0,0,_,_] * (B_p * PXC_hc[:,_,_,:,:,:,:]).sum(axis=(5,6))
    dXi = aggr_inc(inc,P['foi_mode'],axis=(0,3,4))
  elif P['foi_mode'] in ['py']:
    B_p = 1 - ((1 - beta) ** A_ap[:,:,_,_,_,_,_,_]).prod(axis=0)
    QA = mix[:,:,:,:,:,_,_] / X.sum(axis=(2,3,4))[_,:,:,_,_,_,_]
    inc = 1 - ((1 - B_p * PXC_hc[:,_,_,:,:,:,:]) ** QA).prod(axis=(5,6))
    dXi = aggr_inc(inc,P['foi_mode'],axis=(0,3,4),Xsus=X[:,:,0,0,0])
  # all non-base cases
  dX[:,:,0,0,0] -= dXi
  dX[:,:,0,1,0] += dXi
  return inc

#@profile
def aggr_inc(inc,foi_mode,axis,Xsus=None,Xinf=None,keepdims=False):
  # returns absolute infections after appropriately aggregating "inc"
  # Xsus only needed if 'foi_mode' in ['py']
  if foi_mode in ['base','lin','rd','ry']:
    return inc.sum(axis=axis,keepdims=keepdims)
  if foi_mode in ['py']:
    return (1 - (1 - inc).prod(axis=axis,keepdims=keepdims)) * Xsus

#@profile
def fix_XK(X,P):
  # we cannot remove more partnerships than we have: X[:,:,1:] <= XKm
  # if this happens (due to turnover) then move the extra X[:,:,1:] to X[:,:,0]
  XKm = X.sum(axis=2,keepdims=True) * np.moveaxis(P['K_psi'],0,2)[...,_]
  XKe = np.maximum(0, X[:,:,1:] - XKm)
  X[:,:,1:] -= XKe
  X[:,:,0 ] += XKe.sum(axis=2)
