import numpy as np
from utils import _,deco,dtfun

labels = {
  'NX':         'Population Size, Absolute (\'000s)',
  'Psi':        'Population Size, Relative',
  'Ph':         'Proportion of HIV+',
  'prevalence': 'Prevalence',
  'incidence':  'Incidence (per person-year)',
  'cuminfect':  'Cumulative Infections (\'000s)',
  'diagnosed':  'Diagnosed among HIV+',
  'treated_u':  'Treated among HIV+',
  'treated_c':  'Treated among Diagnosed',
  'vls_u':      'Virally Suppressed among HIV+',
  'vls_c':      'Virally Suppressed among Treated',
  'dx_rate':    'Diagnosis rate (per person-year)',
  'tx_rate':    'Treatment rate (per person-year)',
  'condom':     'Condom Use (Proportion of Acts)',
  'circum':     'Proportion of Men Circumcised',
  'gud':        'Relative GUD Prevalence',
}

def by_name(name):
  # get a function in this module by its name
  return globals()[name]

@deco.nanzero
def vs_pop(ofun,R,pop1,pop2,vsop='1/2',aggr=True,**kwds):
  # compute ofun from R for pop1 and pop2, then do some operation to combine the results
  if isinstance(ofun,str): ofun = by_name(ofun)
  O1 = ofun(R,**pop1,**kwds,aggr=aggr)
  O2 = ofun(R,**pop2,**kwds,aggr=aggr)
  if   vsop=='1/2':    return (O1/O2)
  elif vsop=='1-2':    return (O1-O2)
  elif vsop=='1-2/1':  return (O1-O2)/O1
  elif vsop=='1-2/2':  return (O1-O2)/O2
  else: raise NotImplementedError('out.vs_pop(): vsop = '+str(vsop))

@deco.nanzero
def vs_R(ofun,R1,R2,vsop='1-2',aggr=True,**kwds):
  if isinstance(ofun,str): ofun = by_name(ofun)
  O1 = ofun(R1,**kwds,aggr=aggr)
  O2 = ofun(R2,**kwds,aggr=aggr)
  if   vsop=='1/2':    return (O1/O2)
  elif vsop=='1-2':    return (O1-O2)
  elif vsop=='1-2/1':  return (O1-O2)/O1
  elif vsop=='1-2/2':  return (O1-O2)/O2
  else: raise NotImplementedError('out.vs_R(): vsop = '+str(vsop))

def vs_label(lab1,lab2,vsop):
  if   vsop=='1/2':    return '{} / {}'.format(lab1,lab2)
  elif vsop=='1-2':    return '{} - {}'.format(lab1,lab2)
  elif vsop=='1-2/1':  return '{} - {} (Rel)'.format(lab1,lab2)
  elif vsop=='1-2/2':  return '{} - {} (Rel)'.format(lab1,lab2)
  else: raise NotImplementedError('out.vs_label(): vsop = '+str(vsop))

# X dimensions: t,s,i,h,c

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
def NX(X,s=None,i=None,aggr=False):
  X = X.sum(axis=(3,4))
  X = X_by_si(X,s=s,i=i)
  return X.sum(axis=(1,2)) if aggr else X

@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def Psi(X,s=None,i=None,aggr=False):
  X  = X.sum(axis=(3,4))
  XS = X.sum(axis=(1,2),keepdims=True)
  X  = X_by_si(X,s=s,i=i)
  return aggratio(X,XS,aggr)

@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def prevalence(X,s=None,i=None,aggr=False):
  X  = X_by_si(X,s=s,i=i)
  XS = X.sum(axis=(3,4)) # <- denominator; numerator -> (1 - susceptible)
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4))
  return aggratio(Xhiv,XS,aggr)

@deco.rmap(args=['X','inc'])
@deco.tslice(targs=['X','inc'])
def incidence(X,inc,s=None,i=None,aggr=False):
  X = X[:,:,:,0,0] # only susceptible
  Xinc = X_by_si(X*inc,s=s,i=i)
  Xsus = X_by_si(X,s=s,i=i)
  return aggratio(Xinc,Xsus,aggr)

@deco.rmap(args=['X','inc'])
def cuminfect(X,inc,tvec,s=None,i=None,aggr=False):
  X  = X[:,:,:,0,0] # only susceptible
  dt = dtfun(tvec)
  Xinc = X_by_si(X*inc,s=s,i=i)
  return np.cumsum(Xinc.sum(axis=(1,2)) * dt if aggr else Xinc * dt[:,_,_], axis=0)

@deco.nanzero
@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def Ph(X,h,s=None,i=None,aggr=False):
  X    = X_by_si(X,s=s,i=i)
  Xh   = (X[:,:,:,_,h] if isinstance(h,int) else X[:,:,:,h]).sum(axis=(3,4))
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4))
  return aggratio(Xh,Xhiv,aggr)

@deco.nanzero
@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def diagnosed(X,s=None,i=None,aggr=False):
  X = X_by_si(X,s=s,i=i)
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4)) # HIV+
  Xdia = X[:,:,:,1:,1:].sum(axis=(3,4)) # diagnosed
  return aggratio(Xdia,Xhiv,aggr)

@deco.nanzero
@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def treated(X,s=None,i=None,aggr=False,cond=False):
  X = X_by_si(X,s=s,i=i)
  if cond:
    Xref = X[:,:,:,1:,1:].sum(axis=(3,4)) # diagnosed
  else:
    Xref = X[:,:,:,1:,:].sum(axis=(3,4)) # HIV+
  Xtre = X[:,:,:,1:,3:].sum(axis=(3,4))
  return aggratio(Xtre,Xref,aggr)

@deco.nanzero
@deco.rmap(args=['X'])
@deco.tslice(targs=['X'])
def vls(X,s=None,i=None,aggr=False,cond=False):
  X = X_by_si(X,s=s,i=i)
  if cond:
    Xref = X[:,:,:,1:,3:].sum(axis=(3,4)) # treated
  else:
    Xref = X[:,:,:,1:,:].sum(axis=(3,4)) # HIV+
  Xvls = X[:,:,:,1:,4:].sum(axis=(3,4))
  return aggratio(Xvls,Xref,aggr)

treated_u = lambda *a,**k: treated(*a,**k,cond=False)
treated_c = lambda *a,**k: treated(*a,**k,cond=True)
vls_u     = lambda *a,**k: vls(*a,**k,cond=False)
vls_c     = lambda *a,**k: vls(*a,**k,cond=True)

@deco.nanzero
@deco.rmap(args=['X','dx_t'])
@deco.tslice(targs=['X','dx_t'])
def dx_rate(X,dx_t,s=None,i=None,aggr=False):
  X = X[:,:,:,1:,0].sum(axis=(3)) # undiagnosed only
  Xdx = X_by_si(X*np.squeeze(dx_t),s=s,i=i)
  XS  = X_by_si(X,s=s,i=i)
  return aggratio(Xdx,XS,aggr)

@deco.nanzero
@deco.rmap(args=['X','tx','Rtx_ht'])
@deco.tslice(targs=['X','tx','Rtx_ht'])
def tx_rate(X,tx,Rtx_ht,s=None,i=None,aggr=False):
  X = X[:,:,:,1:6,1]
  Xtx = X_by_si(X*tx*Rtx_ht,s=s,i=i).sum(axis=3)
  XS  = X_by_si(X,s=s,i=i).sum(axis=3)
  return aggratio(Xtx,XS,aggr)

@deco.nanzero
@deco.rmap(args=['X','P'])
@deco.tslice(targs=['X'])
def X_rate(X,P,rate,s=None,i=None,aggr=False):
  Xrate = X_by_si(X*P[rate][_,],s=s,i=i).sum(axis=(3,4))
  XS    = X_by_si(X,s=s,i=i).sum(axis=(3,4))
  return aggratio(Xrate,XS,aggr)

@deco.rmap(args=['PA_condom_t'])
@deco.tslice(targs=['PA_condom_t'])
def condom(PA_condom_t,p,aggr=None):
  return np.squeeze(PA_condom_t)[:,p]

@deco.rmap(args=['PA_circum_t'])
@deco.tslice(targs=['PA_circum_t'])
def circum(PA_circum_t,aggr=None):
  return np.squeeze(PA_circum_t)[:]

@deco.rmap(args=['P_gud_t'])
@deco.tslice(targs=['P_gud_t'])
def gud(P_gud_t,aggr=None):
  return np.squeeze(P_gud_t)[:]
