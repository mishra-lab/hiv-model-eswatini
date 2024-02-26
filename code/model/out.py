# functions for computing (stratified) model outputs

import numpy as np
from itertools import product as iprod
from utils import _,deco,xdi,dtfun,itslice
from model import system,foi,strats

# ------------------------------------------------------------------------------
# helper data & functions

labels = {
  'NX':         'Population size, absolute (\'000s)',
  'Psi':        'Population size, relative',
  'Ph':         'Proportion of PLHIV',
  'prevalence': 'Prevalence',
  'incidence':  'Incidence (per person-year)',
  'cuminfect':  'Cumulative infections (\'000s)',
  'tdsc':       'Transmission-driven seroconcordance',
  'diagnosed':  'Diagnosed among PLHIV',
  'treated_u':  'Treated among PLHIV',
  'treated_c':  'Treated among diagnosed',
  'vls_u':      'VLS among PLHIV',
  'vls_c':      'VLS among treated',
  'dx_rate':    'Diagnosis rate (per person-year)',
  'tx_rate':    'Treatment rate (per person-year)',
  'condom':     'Condom use (proportion of acts)',
  'circum':     'Proportion of men circumcised',
  'gud':        'Relative GUD prevalence',
}

qs = [0,.025,.05,.1,.25,.4,.45,.475,.5,.525,.55,.6,.75,.9,.95,.975,1]

def by_name(name):
  # get a function in this module by its name
  return globals()[name]

@deco.nanzero
def vs_fun(O1,O2,vsop):
  if   vsop=='1/2':    return (O1/O2)
  elif vsop=='2/1':    return (O2/O1)
  elif vsop=='1-2':    return (O1-O2)
  elif vsop=='2-1':    return (O2-O1)
  elif vsop=='1-2/1':  return (O1-O2)/O1
  elif vsop=='2-1/1':  return (O2-O1)/O1
  elif vsop=='1-2/2':  return (O1-O2)/O2
  elif vsop=='2-1/2':  return (O2-O1)/O2
  else: raise NotImplementedError('out.vs_fun(): vsop = '+str(vsop))

def vs_label(lab1,lab2,vsop):
  if   vsop=='1/2':    return '{} / {}'.format(lab1,lab2)
  elif vsop=='2/1':    return '{} / {}'.format(lab2,lab1)
  elif vsop=='1-2':    return '{} - {}'.format(lab1,lab2)
  elif vsop=='2-1':    return '{} - {}'.format(lab2,lab1)
  elif vsop=='1-2/1':  return '{} - {} (REF=1)'.format(lab1,lab2)
  elif vsop=='2-1/1':  return '{} - {} (REF=1)'.format(lab2,lab1)
  elif vsop=='1-2/2':  return '{} - {} (REF=2)'.format(lab1,lab2)
  elif vsop=='2-1/2':  return '{} - {} (REF=2)'.format(lab2,lab1)
  if   vsop=='\n':     return '{}\nvs. {}'.format(lab1,lab2) # HACK
  else: raise NotImplementedError('out.vs_label(): vsop = '+str(vsop))

@deco.nanzero
def vs_ind(ofun,R,ind1,ind2,vsop='1/2',aggr=True,**kwds):
  # compute ofun from R for ind1 & ind2, then combine via vs_fun
  if isinstance(ofun,str): ofun = by_name(ofun)
  O1 = ofun(R,**ind1,**kwds,aggr=aggr)
  O2 = ofun(R,**ind2,**kwds,aggr=aggr)
  return vs_fun(O1,O2,vsop)

@deco.nanzero
def vs_R(ofun,R1,R2,vsop='1-2',aggr=True,**kwds):
  # compute ofun for R1 & R2 with same kwds, then combine via vs_fun
  if isinstance(ofun,str): ofun = by_name(ofun)
  O1 = ofun(R1,**kwds,aggr=aggr)
  O2 = ofun(R2,**kwds,aggr=aggr)
  return vs_fun(O1,O2,vsop)

def aggratio(X1,X2,aggr,axis=None):
  # helper function: many outputs are defined as X1 / X2
  # sometimes we need the elementwise result
  # sometimes we need X1.sum() / X2.sum()
  if axis is None: axis = (1,2)
  if aggr:
    return X1.sum(axis=axis) / X2.sum(axis=axis)
  else:
    return X1 / X2

# X.shape: (t:*, s:2, i:4, h:6, c:5)
# Xk.shape: (t:*, s:2, i:4, k:5, h:6, c:5) - tdsc only

# ------------------------------------------------------------------------------
# main output functions

@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def NX(X,s=None,i=None,aggr=True):
  X = X.sum(axis=(3,4)) # sum health & care (all)
  X = xdi(X,{1:s,2:i})  # select sex & activity
  return X.sum(axis=(1,2)) if aggr else X

@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def Psi(X,s=None,i=None,aggr=True):
  X  = X.sum(axis=(3,4)) # sum health & care (all)
  XS = xdi(X,{1:_,2:_})  # sum sex & activity (denom)
  X  = xdi(X,{1:s,2:i})  # select sex & activity (num)
  return aggratio(X,XS,aggr)

@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def prevalence(X,s=None,i=None,aggr=True):
  X  = xdi(X,{1:s,2:i})  # select sex & activity
  XS = X.sum(axis=(3,4)) # sum health & care (denom)
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4)) # select plhiv, sum health & care (num)
  return aggratio(Xhiv,XS,aggr)

@deco.rmap(Rk=['X','inc'],Pk=['foi_mode'])
@deco.tslice(tk=['X','inc'])
def incidence(X,inc,foi_mode,s=None,i=None,aggr=True):
  # total infections (per year); uses foi.aggr_inc due to FOI cases
  inf = foi.aggr_inc(inc,foi_mode,axis=(1,4,5),Xsus=X[:,:,:,0,0])
  sus_si = xdi(X[:,:,:,0,0],{1:s,2:i}) # select sex & activity among sus (denom)
  inf_si = xdi(inf,{1:s,2:i})          # select sex & activity new infs (num)
  return aggratio(inf_si,sus_si,aggr)

@deco.rmap(Rk=['X','inc'],Pk=['foi_mode'])
def cuminfect(X,inc,foi_mode,tvec,s=None,i=None,aggr=True,t0=None):
  dt = dtfun(tvec) # timestep sizes
  # total infections (per year); uses foi.aggr_inc due to FOI cases
  inf = foi.aggr_inc(inc,foi_mode,axis=(1,4,5),Xsus=X[:,:,:,0,0])
  inf_si = xdi(inf,{1:s,2:i}) # select sex & activity new infs
  # sum new infs across sex & activity if aggr, mult by dt
  inf_dt = inf_si.sum(axis=(1,2)) * dt if aggr else inf_si * dt[:,_,_]
  if t0: # zero new infs before t0
    inf_dt[tvec < t0] = 0
  return np.cumsum(inf_dt,axis=0)

@deco.nanzero
@deco.rmap(Rk=['Xk'],Pk=['K_psi'])
@deco.tslice(tk=['Xk'])
def tdsc(Xk,K_psi,p=None,s=None,i=None,aggr=True,sus=False):
  # tdsc = transmission-driven seroconcordance
  # only works for foi_mode='base' & system.run(...,Xk=True)
  # Xk.shape: (t:*, s:2, i:4, k:5, h:6, c:5)
  # select inf+sus or inf only, reshape, sum health & care
  Xk = np.moveaxis(Xk if sus else Xk[:,:,:,:,1:,:],3,1).sum(axis=(4,5))
  # sum seroconc dim, mult by # ptrs, select ptr-type, sex, & activity (denom)
  XK   = xdi(xdi(Xk,{1:_}) * np.squeeze(K_psi)[_,:,:,:],{1:p,2:s,3:i})
  # select seroconc only, ptr-type, sex, & activity (num)
  XKsc = xdi(Xk[:,1:,:,:],{1:p,2:s,3:i})
  return aggratio(XKsc,XK,aggr,axis=(1,2,3))

def can_tdsc(R):
  # check if we can compute tdsc from this R
  return ('Xk' in R) and (R['P']['foi_mode'] == 'base')

@deco.nanzero
@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def Ph(X,h,s=None,i=None,aggr=True):
  X    = xdi(X,{1:s,2:i})              # select sex & activity (all)
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4)) # select plhiv, sum health & care (denom)
  Xh   = xdi(X,{3:h}).sum(axis=(3,4))  # select health, sum health & care (num)
  return aggratio(Xh,Xhiv,aggr)

@deco.nanzero
@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def diagnosed(X,s=None,i=None,aggr=True):
  X = xdi(X,{1:s,2:i})                  # select sex & activity (all)
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4))  # select plhiv, sum health & care (denom)
  Xdia = X[:,:,:,1:,1:].sum(axis=(3,4)) # select diagnosed, sum health & care (num)
  return aggratio(Xdia,Xhiv,aggr)

@deco.nanzero
@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def treated(X,s=None,i=None,aggr=True,cond=False):
  X = xdi(X,{1:s,2:i})                    # select sex & activity (all)
  if cond:
    Xref = X[:,:,:,1:,1:].sum(axis=(3,4)) # select diagnosed, sum health & care (denom)
  else:
    Xref = X[:,:,:,1:,:].sum(axis=(3,4))  # select plhiv, sum health & care (denom)
  Xtre = X[:,:,:,1:,3:].sum(axis=(3,4))   # select treated, sum health & care (num)
  return aggratio(Xtre,Xref,aggr)

@deco.nanzero
@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def vls(X,s=None,i=None,aggr=True,cond=False):
  X = xdi(X,{1:s,2:i})                    # select sex & activity (all)
  if cond:
    Xref = X[:,:,:,1:,3:].sum(axis=(3,4)) # select treated, sum health & care (denom)
  else:
    Xref = X[:,:,:,1:,:].sum(axis=(3,4))  # select plhiv, sum health & care (denom)
  Xvls = X[:,:,:,1:,4:].sum(axis=(3,4))   # select vls, sum health & care (num)
  return aggratio(Xvls,Xref,aggr)

treated_u = lambda *a,**k: treated(*a,**k,cond=False)
treated_c = lambda *a,**k: treated(*a,**k,cond=True)
vls_u     = lambda *a,**k: vls(*a,**k,cond=False)
vls_c     = lambda *a,**k: vls(*a,**k,cond=True)

@deco.nanzero
@deco.rmap(Rk=['X','dx_sit'])
@deco.tslice(tk=['X','dx_sit'])
def dx_rate(X,dx_sit,s=None,i=None,aggr=True):
  X = X[:,:,:,1:,0].sum(axis=(3))           # select undiagnosed, sum health
  XS  = xdi(X,{1:s,2:i})                    # select sex & activity (denom)
  Xdx = xdi(X*np.squeeze(dx_sit),{1:s,2:i}) # mult by rate, select sex & activity (num)
  return aggratio(Xdx,XS,aggr)

@deco.nanzero
@deco.rmap(Rk=['X','tx_sit','Rtx_ht'])
@deco.tslice(tk=['X','tx_sit','Rtx_ht'])
def tx_rate(X,tx_sit,Rtx_ht,s=None,i=None,aggr=True):
  X = X[:,:,:,1:,1]                   # select diagnosed
  XS  = xdi(X,{1:s,2:i}).sum(axis=3)  # select sex & activity, sum health (denom)
  # mult by rate, select sex & activity, sum health (num)
  Xtx = xdi(X*np.squeeze(tx_sit*Rtx_ht),{1:s,2:i}).sum(axis=3)
  return aggratio(Xtx,XS,aggr)

@deco.rmap(Rk=['PF_condom_t'])
@deco.tslice(tk=['PF_condom_t'])
def condom(PF_condom_t,p,aggr=None):
  return np.squeeze(PF_condom_t)[:,p] # select ptr-type

@deco.rmap(Rk=['PF_circum_t'])
@deco.tslice(tk=['PF_circum_t'])
def circum(PF_circum_t,aggr=None):
  return np.squeeze(PF_circum_t)[:]

@deco.rmap(Rk=['X','inc'],Pk=['foi_mode'])
@deco.tslice(tk=['X','inc'])
def infections(X,inc,foi_mode,p,fs,fi,ts,ti):
  # total infections (per year); uses foi.aggr_inc due to FOI cases
  # NOTE: p,fs,fi,ts,ti must be single values!
  return foi.aggr_inc(inc[:,p,ts,ti,fs,fi],foi_mode,axis=(),Xsus=X[:,ts,ti,0,0])

# ------------------------------------------------------------------------------
# output collection functions

def wiw(R1s,tvec,t,R2s=None,vsop='1-2'):
  # who infected whom: tabular data (list of lists): infections per year
  # stratified by time/ptr-type/to-group/from-group
  # optionally get (paired) difference between R1s & R2s
  # use quantiles (vs for every R/id) to reduce data size
  if isinstance(R1s,dict): R1s = [R1s]
  if isinstance(R2s,dict): R2s = [R2s]
  aggrop = lambda inf: np.nanquantile(inf,qs,axis=0)
  cols = ['q'+str(q) for q in qs]
  data = [['t','p','fs','fi','ts','ti']+cols] # init cols
  kwds = dict(tvec=tvec,t=t)
  grid = dict(p=range(4),fs=range(2),fi=range(4),ts=range(2),ti=range(4))
  # t=time, p=ptr-type, f*=from, t*=to, *s=sex, *i=activity
  for p,fs,fi,ts,ti in iprod(*grid.values()):
    kwds.update(p=p,fs=fs,fi=fi,ts=ts,ti=ti)
    inf = aggrop([infections(R1,**kwds) for R1 in R1s]) if R2s is None else \
          aggrop([vs_fun(infections(R1,**kwds),infections(R2,**kwds),vsop) for R1,R2 in zip(R1s,R2s)])
    # append rows (all time points) for this combination of p,fs,fi,ts,ti
    data += [[tk,p,fs,fi,ts,ti]+inf[:,k].tolist() for k,tk in enumerate(t)]
  return data

def expo(R1s,tvec,t,onames,skeys,R2s=None,vsop='raw',ecols=None,mode='q',**kwds):
  # tabular data (dict of lists): any outputs, for same strata
  # optionally get (paired) difference between R1s & R2s
  if mode == 'q': # quantiles
    aggrop = lambda os: np.nanquantile(os,qs,axis=0)
    cols = ['q'+str(q) for q in qs]
  if mode == 'id': # every R/id
    aggrop = lambda os: np.array(os)
    cols = ['i'+str(R['P']['id']) for R in R1s]
  sg,og,tg = [g.flatten().tolist() for g in np.meshgrid(skeys,onames,t)] # grid
  if ecols is None: ecols = {} # option to add extra columns (static values)
  ecols.update(op=vsop)
  E = dict(out=og,pop=sg,t=tg, # init cols: grid
    **{k:[v]*len(tg) for k,v in ecols.items()}, # extra cols
    **{k:[] for k in cols}) # empty data cols
  for oname in onames: # outputs
    ofun = by_name(oname)
    for skey in skeys: # strata
      if oname == 'cuminfect': # special case: cannot use @deco.tslice
        sofun = lambda R: ofun(R,**strats[skey].ind,tvec=tvec,**kwds)[itslice(t,tvec)]
      else:
        sofun = lambda R: ofun(R,**strats[skey].ind,tvec=tvec,t=t,**kwds)
      # compute the output
      if R2s is None:
        osx = aggrop([sofun(R) for R in R1s])
      else:
        osx = aggrop([vs_fun(sofun(R1),sofun(R2),vsop) for R1,R2 in zip(R1s,R2s)])
      # TODO: below might be slow, vectorize with np?
      for i,col in enumerate(cols): # append to data columns
        E[col] += osx[i,:].tolist()
  return E
