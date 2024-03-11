# functions for computing (stratified) model outputs

import numpy as np
from utils import _,deco,linear_comb
from model import tol

foi_modes = [
  'lin',  # linear          {unused}
  'rd',   # rate duration   {partnership-duration}
  'ry',   # rate year       {partnership-year}
  'py',   # proportion year {all partnerships per year}
  'base', # EPA (default)   {proposed approach}
]

#@profile
def get_beta(P,t):
  # beta: probability of transmission per sex act
  # beta.shape = (a:2, p:4, s:2, i:4, s':2, i':4, h':6, c':5)
  RbA_condom = linear_comb(P['PF_condom_t'](t) * P['RPF_condom_a'], P['Rbeta_condom'], 1)
  RbA_circum = linear_comb(P['PF_circum_t'](t), P['Rbeta_circum'], 1)
  P_gud_t = P['P_gud'] * P['RP_gud_t'](t) # GUD = genital ulcer diseases
  Rbeta_gud_sus = linear_comb(P_gud_t,1+P['aRbeta_gud_sus'],1).reshape([1,1,2,4,1,1,1,1]) # self
  Rbeta_gud_inf = linear_comb(P_gud_t,1+P['aRbeta_gud_inf'],1).reshape([1,1,1,1,2,4,1,1]) # other
  # product of base prob & relative factors, force beta <= 0.5
  return np.minimum(.5, P['beta_a'] * Rbeta_gud_sus * Rbeta_gud_inf * RbA_condom * RbA_circum)

@deco.nowarn
#@profile
def get_mix(XC,P):
  # mix: population-scale mixing (total partners), mix.shape = (p:4, s:2, i:4, s':2, i':4)
  # XC: total partners, XC.shape = (p:4, s:2, i:4)
  # M0: random mixing, M0.shape = (p:4, i:4 {women}, i':4 {men})
  M0 = XC[:,0,:,_] * XC[:,1,_,:] / XC.sum(axis=2).mean(axis=1)[:,_,_] + tol/10
  m1 = M0.sum(axis=1) # total for men
  m2 = M0.sum(axis=2) # total for women
  # print(m1 / XC[:,1,:]) # DEBUG == 1, unless XC unbalanced
  # print(m2 / XC[:,0,:]) # DEBUG == 1, unless XC unbalanced
  M = M0 * np.exp(P['pref_pii']) # apply mixing log-odds
  # iterative proportional fitting
  for k in range(100):
    r1 = m1 / M.sum(axis=1)
    M *= r1[:,_,:]
    r2 = m2 / M.sum(axis=2)
    M *= r2[:,:,_]
    if (abs(r1-1) < tol).all() and (abs(r2-1) < tol).all():
      break # close enough, usually k < 20
  M[abs(M)<tol] = 0 # fix rounding errors
  P['mix'][:,0,:,1,:] = M
  P['mix'][:,1,:,0,:] = M.swapaxes(1,2)
  return P['mix']

@deco.nowarn
#@profile
def get_apply_inc(dX,X,t,P):
  # inc.shape = (p:4, s:2, i:4, s':2, i':4)
  # > if foi_mode in ['base','lin','rd','ry']: return *absolute* infections (not per susceptible)
  # > if foi_mode in ['py']: return *probability* of infection (aggr must be deferred)
  # C_psik = partner numbers (K) or rates of ptr change (Q); A = sex acts per partnership
  # C_psik.shape = (p:4, s:2, i:4, k:5)
  if P['foi_mode'] in ['base']:
    C_psik = P['K_psi'] - P['aK_pk'] # K (ptr count), EPA adjusted
    fix_XK(X,P) # see function comments
  elif P['foi_mode'] in ['lin']:
    C_psik = P['K_psi'] # K (ptr count), no adjustment
  elif P['foi_mode'] in ['rd']:
    C_psik = P['K_psi'] / P['dur_p'][:,_,_,_] # Q (ptr rate)
    A_ap = P['F_ap'] * P['dur_p'] # A (sex acts, full duration)
  elif P['foi_mode'] in ['ry','py']:
    C_psik = P['K_psi'] / P['dur_p_1'][:,_,_,_] # Q_1 (ptr rate >= 1)
    A_ap = P['F_ap'] * P['dur_p_1'] # A (sex acts, duration <= 1 year)
  # setup mixing & compute prevalence
  XC_psihc = (X[_,:,:,:,:,:] * C_psik[:,:,:,:,_,_]).sum(axis=3) # total effective partners
  XC_psihc[XC_psihc<0] = 0 # fix rounding errors
  XC_psi = XC_psihc.sum(axis=(3,4)) # shape = (p:4, s:2, i:4)
  # Phc = % in strata (h,c) among (p,s,i); unless foi_mode = 'base', Phc_XC_psi = Phc_X_psi
  Phc_XC_psi = XC_psihc / (XC_psi[:,:,:,_,_] + tol/10) # shape = (p:4, s:2, i:4, h:6, c:5)
  # compute population-scale mixing
  mix = get_mix(XC_psi,P) * P['mix_mask'] # shape = (p:4, s:2, i:4, s':2, i':4)
  # compute per-act probability
  beta = get_beta(P,t) # shape = (a:2, p:4, s:2, i:4, s':2, i':4, h':6, c':5)
  # compute & apply force of infection
  if P['foi_mode'] in ['base']:
    # inc = {# ptrs} * {% sus} * {% inf} * {beta per-act} * {act freq}
    Fbeta = (beta * P['F_ap'][:,:,_,_,_,_,_,_]).sum(axis=0) # sum sex act types (a:2)
    # inc.shape = (p:4, s:2, i:4, s':2, i':4, h':6, c':5)
    inc = mix[:,:,:,:,:,_,_] * Phc_XC_psi[:,:,:,0,0,_,_,_,_] * Fbeta * Phc_XC_psi[:,_,_,:,:,:,:]
    dXi = inc.sum(axis=(3,4,5,6)) # acquisition: (p:4, s:2, i:4)
    dX[:,:,0 ,0,0] -= dXi.sum(axis=0)
    dX[:,:,1:,1,0] += np.moveaxis(dXi,0,2)
    dXi = inc.sum(axis=(1,2)) # transmission: (p:4, s':2, i':4, h':6, c':5)
    dX[:,:,0 ,:,:] -= dXi.sum(axis=0)
    dX[:,:,1:,:,:] += np.moveaxis(dXi,0,2)
    dXi = X[:,:,1:,:,:] / P['dur_p'][_,_,:,_,_] # new ptrs: (s:2, i:4, k:4, h:6, c:5)
    dX[:,:,1:,:,:] -= dXi
    dX[:,:,0 ,:,:] += dXi.sum(axis=2)
    return inc.sum(axis=(5,6)) # done
  elif P['foi_mode'] in ['lin']:
    # inc = {# ptrs} * {% sus} * {% inf} * {beta per-act} * {act freq}
    Fbeta = (beta * P['F_ap'][:,:,_,_,_,_,_,_]).sum(axis=0) # sum acts
    inc = mix * Phc_XC_psi[:,:,:,0,0,_,_] * (Fbeta * Phc_XC_psi[:,_,_,:,:,:,:]).sum(axis=(5,6))
    dXi = aggr_inc(inc,P['foi_mode'],axis=(0,3,4)) # sum across ptrs
  elif P['foi_mode'] in ['rd','ry']:
    # inc = {# ptrs} * {% sus} * {% inf} * (1 - (1 - {beta per-act}) ^ {acts per-ptr})
    B_p = 1 - ((1 - beta) ** A_ap[:,:,_,_,_,_,_,_]).prod(axis=0) # prod acts
    inc = mix * Phc_XC_psi[:,:,:,0,0,_,_] * (B_p * Phc_XC_psi[:,_,_,:,:,:,:]).sum(axis=(5,6))
    dXi = aggr_inc(inc,P['foi_mode'],axis=(0,3,4)) # sum across ptrs
  elif P['foi_mode'] in ['py']:
    # B = (1 - (1 - {beta per-act}) ^ {acts per-ptr})
    # inc = {% sus} * (1 - (1 - {B} * {% inf}) ^ {# ptrs})
    B_p = 1 - ((1 - beta) ** A_ap[:,:,_,_,_,_,_,_]).prod(axis=0) # prod acts
    mix_pp = mix[:,:,:,:,:,_,_] / X.sum(axis=(2,3,4))[_,:,:,_,_,_,_]
    inc = 1 - ((1 - B_p * Phc_XC_psi[:,_,_,:,:,:,:]) ** mix_pp).prod(axis=(5,6))
    dXi = aggr_inc(inc,P['foi_mode'],axis=(0,3,4),Xsus=X[:,:,0,0,0]) # prod across ptrs
  # all non-base cases
  dX[:,:,0,0,0] -= dXi # sus
  dX[:,:,0,1,0] += dXi # inf: acute & undx
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
