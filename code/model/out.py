import numpy as np
from itertools import product as iprod
from utils import _,deco,dtfun,itslice
from model import system,foi,slicers
# TODO: integrate slicers fully?

labels = {
  'NX':         'Population size, absolute (\'000s)',
  'Psi':        'Population size, relative',
  'Ph':         'Proportion of PLHIV',
  'prevalence': 'Prevalence',
  'incidence':  'Incidence (per person-year)',
  'cuminfect':  'Cumulative infections (\'000s)',
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
  else: raise NotImplementedError('out.vs_pop(): vsop = '+str(vsop))

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
def vs_pop(ofun,R,pop1,pop2,vsop='1/2',aggr=True,**kwds):
  # compute ofun from R for pop1 and pop2, then do some operation to combine the results
  if isinstance(ofun,str): ofun = by_name(ofun)
  O1 = ofun(R,**pop1,**kwds,aggr=aggr)
  O2 = ofun(R,**pop2,**kwds,aggr=aggr)
  return vs_fun(O1,O2,vsop)

@deco.nanzero
def vs_R(ofun,R1,R2,vsop='1-2',aggr=True,**kwds):
  if isinstance(ofun,str): ofun = by_name(ofun)
  O1 = ofun(R1,**kwds,aggr=aggr)
  O2 = ofun(R2,**kwds,aggr=aggr)
  return vs_fun(O1,O2,vsop)

# X.shape: (t:*, s:2, i:4, h:6, c:5)

def aggratio(X1,X2,aggr,axis=None):
  if axis is None: axis = (1,2)
  if aggr:
    return X1.sum(axis=axis) / X2.sum(axis=axis)
  else:
    return X1 / X2

def X_by_si(X,s=None,i=None):
  X = X.sum(axis=1,keepdims=True) if s is None \
      else X[:,_,s] if isinstance(s,int) \
      else X[:,s]
  X = X.sum(axis=2,keepdims=True) if i is None \
      else X[:,:,_,i] if isinstance(i,int) \
      else X[:,:,i]
  return X

@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def NX(X,s=None,i=None,aggr=True):
  X = X.sum(axis=(3,4))
  X = X_by_si(X,s=s,i=i)
  return X.sum(axis=(1,2)) if aggr else X

@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def Psi(X,s=None,i=None,aggr=True):
  X  = X.sum(axis=(3,4))
  XS = X.sum(axis=(1,2),keepdims=True)
  X  = X_by_si(X,s=s,i=i)
  return aggratio(X,XS,aggr)

@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def prevalence(X,s=None,i=None,aggr=True):
  X  = X_by_si(X,s=s,i=i)
  XS = X.sum(axis=(3,4)) # <- denominator; numerator -> (1 - susceptible)
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4))
  return aggratio(Xhiv,XS,aggr)

@deco.rmap(args=['X','inc','foi_mode'])
@deco.tslice(targs=['X','inc'])
def incidence(X,inc,foi_mode,s=None,i=None,aggr=True):
  inf = foi.aggr_inc(inc,foi_mode,axis=(1,4,5),Xsus=X[:,:,:,0,0])
  inf_si = X_by_si(inf,s=s,i=i)
  sus_si = X_by_si(X[:,:,:,0,0],s=s,i=i)
  return aggratio(inf_si,sus_si,aggr)

@deco.rmap(args=['X','inc','foi_mode'])
def cuminfect(X,inc,foi_mode,tvec,s=None,i=None,aggr=True,t0=None):
  dt = dtfun(tvec)
  inf = foi.aggr_inc(inc,foi_mode,axis=(1,4,5),Xsus=X[:,:,:,0,0])
  inf_si = X_by_si(inf,s=s,i=i)
  inf_dt = inf_si.sum(axis=(1,2)) * dt if aggr else inf_si * dt[:,_,_]
  if t0:
    inf_dt[tvec < t0] = 0
  return np.cumsum(inf_dt,axis=0)

@deco.nanzero
@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def Ph(X,h,s=None,i=None,aggr=True):
  X    = X_by_si(X,s=s,i=i)
  Xh   = (X[:,:,:,_,h] if isinstance(h,int) else X[:,:,:,h]).sum(axis=(3,4))
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4))
  return aggratio(Xh,Xhiv,aggr)

@deco.nanzero
@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def diagnosed(X,s=None,i=None,aggr=True):
  X = X_by_si(X,s=s,i=i)
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4)) # PLHIV
  Xdia = X[:,:,:,1:,1:].sum(axis=(3,4)) # diagnosed
  return aggratio(Xdia,Xhiv,aggr)

@deco.nanzero
@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def treated(X,s=None,i=None,aggr=True,cond=False):
  X = X_by_si(X,s=s,i=i)
  if cond:
    Xref = X[:,:,:,1:,1:].sum(axis=(3,4)) # diagnosed
  else:
    Xref = X[:,:,:,1:,:].sum(axis=(3,4)) # PLHIV
  Xtre = X[:,:,:,1:,3:].sum(axis=(3,4))
  return aggratio(Xtre,Xref,aggr)

@deco.nanzero
@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def vls(X,s=None,i=None,aggr=True,cond=False):
  X = X_by_si(X,s=s,i=i)
  if cond:
    Xref = X[:,:,:,1:,3:].sum(axis=(3,4)) # treated
  else:
    Xref = X[:,:,:,1:,:].sum(axis=(3,4)) # PLHIV
  Xvls = X[:,:,:,1:,4:].sum(axis=(3,4))
  return aggratio(Xvls,Xref,aggr)

treated_u = lambda *a,**k: treated(*a,**k,cond=False)
treated_c = lambda *a,**k: treated(*a,**k,cond=True)
vls_u     = lambda *a,**k: vls(*a,**k,cond=False)
vls_c     = lambda *a,**k: vls(*a,**k,cond=True)

@deco.nanzero
@deco.rmap(args=['X','dx_sit'])
@deco.tslice(targs=['X','dx_sit'])
def dx_rate(X,dx_sit,s=None,i=None,aggr=True):
  X = X[:,:,:,1:,0].sum(axis=(3)) # undiagnosed only
  Xdx = X_by_si(X*np.squeeze(dx_sit),s=s,i=i)
  XS  = X_by_si(X,s=s,i=i)
  return aggratio(Xdx,XS,aggr)

@deco.nanzero
@deco.rmap(args=['X','tx_sit','Rtx_ht'])
@deco.tslice(targs=['X','tx_sit','Rtx_ht'])
def tx_rate(X,tx_sit,Rtx_ht,s=None,i=None,aggr=True):
  X = X[:,:,:,1:,1] # diagnosed only
  Xtx = X_by_si(X*np.squeeze(tx_sit*Rtx_ht),s=s,i=i).sum(axis=3)
  XS  = X_by_si(X,s=s,i=i).sum(axis=3)
  return aggratio(Xtx,XS,aggr)

@deco.rmap(args=['PF_condom_t'])
@deco.tslice(targs=['PF_condom_t'])
def condom(PF_condom_t,p,aggr=None):
  return np.squeeze(PF_condom_t)[:,p]

@deco.rmap(args=['PF_circum_t'])
@deco.tslice(targs=['PF_circum_t'])
def circum(PF_circum_t,aggr=None):
  return np.squeeze(PF_circum_t)[:]

@deco.rmap(args=['X','inc','foi_mode'])
@deco.tslice(targs=['X','inc'])
def infections(X,inc,foi_mode,p,fs,fi,ts,ti):
  # NOTE: p,fs,fi,ts,ti must be single values!
  return foi.aggr_inc(inc[:,p,ts,ti,fs,fi],foi_mode,axis=(),Xsus=X[:,ts,ti,0,0])

def wiw(R1s,tvec,t,R2s=None,vsop='1-2'):
  if isinstance(R1s,dict): R1s = [R1s]
  if isinstance(R2s,dict): R2s = [R2s]
  qs = [0,.025,.05,.1,.25,.4,.45,.475,.5,.525,.55,.6,.75,.9,.95,.975,1]
  aggrop = lambda inf: np.nanquantile(inf,qs,axis=0)
  cols = ['q'+str(q) for q in qs]
  data = [['t','p','fs','fi','ts','ti']+cols]
  kwds = dict(tvec=tvec,t=t)
  grid = dict(p=range(4),fs=range(2),fi=range(4),ts=range(2),ti=range(4))
  for p,fs,fi,ts,ti in iprod(*grid.values()):
    kwds.update(p=p,fs=fs,fi=fi,ts=ts,ti=ti)
    inf = aggrop([infections(R1,**kwds) for R1 in R1s]) if R2s is None else \
          aggrop([vs_fun(infections(R1,**kwds),infections(R2,**kwds),vsop) for R1,R2 in zip(R1s,R2s)])
    data += [[tk,p,fs,fi,ts,ti]+inf[:,k].tolist() for k,tk in enumerate(t)]
  return data
