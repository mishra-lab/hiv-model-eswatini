import numpy as np
from utils import _,deco,dtfun
from model import foi # TODO

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

@deco.rmap(args=['X','inc'])
@deco.tslice(targs=['X','inc'])
def whoinfectwhom(X,inc,p=None,fpop=None,tpop=None,aggr=True):
  # total infections between pop1 & pop2 along partnership type p
  # TODO: implement cumulative?
  fs = None if fpop is None else fpop.get('s')
  fi = None if fpop is None else fpop.get('i')
  ts = None if tpop is None else tpop.get('s')
  ti = None if tpop is None else tpop.get('i')
  # aggregating non-self dimensions: (p[1], s'[4], i'[5])
  inc = foi.aggr_inc(inc,axis=1,keepdims=True) if p is None \
        else inc[:,_,p] if isinstance(p,int) else inc[:,p]
  inc = foi.aggr_inc(inc,axis=4,keepdims=True) if fs is None \
        else inc[:,:,:,:,_,fs] if isinstance(fs,int) else inc[:,:,:,:,fs]
  inc = foi.aggr_inc(inc,axis=5,keepdims=True) if fi is None \
        else inc[:,:,:,:,:,_,fi] if isinstance(fi,int) else inc[:,:,:,:,:,fi]
  inf = foi.aggr_inc(inc,axis=(1,4,5),keepdims=True,Xsus=X[:,_,:,:,0,0,_,_]) if aggr \
        else foi.aggr_inc(inc,axis=(),keepdims=True,Xsus=X[:,_,:,:,0,0,_,_])
  inf = inf.sum(axis=2,keepdims=True) if ts is None \
        else inf[:,:,_,ts] if isinstance(ts,int) else inf[:,:,ts]
  inf = inf.sum(axis=3,keepdims=True) if ti is None \
        else inf[:,:,:,_,ti] if isinstance(ti,int) else inf[:,:,:,ti]
  inf = inf.sum(axis=(2,3),keepdims=True) if aggr else inf
  return np.squeeze(inf) if aggr else inf

@deco.rmap(args=['X','inc'])
@deco.tslice(targs=['X','inc'])
def incidence(X,inc,s=None,i=None,aggr=True):
  inf = foi.aggr_inc(inc,axis=(1,4,5),Xsus=X[:,:,:,0,0])
  inf_si = X_by_si(inf,s=s,i=i)
  sus_si = X_by_si(X[:,:,:,0,0],s=s,i=i)
  return aggratio(inf_si,sus_si,aggr)

@deco.rmap(args=['X','inc'])
def cuminfect(X,inc,tvec,s=None,i=None,aggr=True,t0=None):
  dt = dtfun(tvec)
  inf = foi.aggr_inc(inc,axis=(1,4,5),Xsus=X[:,:,:,0,0])
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
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4)) # HIV+
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
    Xref = X[:,:,:,1:,:].sum(axis=(3,4)) # HIV+
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
def dx_rate(X,dx_t,s=None,i=None,aggr=True):
  X = X[:,:,:,1:,0].sum(axis=(3)) # undiagnosed only
  Xdx = X_by_si(X*np.squeeze(dx_t),s=s,i=i)
  XS  = X_by_si(X,s=s,i=i)
  return aggratio(Xdx,XS,aggr)

@deco.nanzero
@deco.rmap(args=['X','tx','Rtx_ht'])
@deco.tslice(targs=['X','tx','Rtx_ht'])
def tx_rate(X,tx,Rtx_ht,s=None,i=None,aggr=True):
  X = X[:,:,:,1:6,1]
  Xtx = X_by_si(X*tx*Rtx_ht,s=s,i=i).sum(axis=3)
  XS  = X_by_si(X,s=s,i=i).sum(axis=3)
  return aggratio(Xtx,XS,aggr)

@deco.nanzero
@deco.rmap(args=['X','P'])
@deco.tslice(targs=['X'])
def X_rate(X,P,rate,s=None,i=None,aggr=True):
  Prate = P[rate][:,:,0,:,:] if P[rate].ndim == 5 else P[rate] # kinda dangerous
  Xrate = X_by_si(X*Prate[_,],s=s,i=i).sum(axis=(3,4))
  XS    = X_by_si(X,s=s,i=i).sum(axis=(3,4))
  return aggratio(Xrate,XS,aggr)

@deco.rmap(args=['PF_condom_t'])
@deco.tslice(targs=['PF_condom_t'])
def condom(PF_condom_t,p,aggr=None):
  return np.squeeze(PF_condom_t)[:,p]

@deco.rmap(args=['PF_circum_t'])
@deco.tslice(targs=['PF_circum_t'])
def circum(PF_circum_t,aggr=None):
  return np.squeeze(PF_circum_t)[:]

@deco.rmap(args=['P_gud_t'])
@deco.tslice(targs=['P_gud_t'])
def gud(P_gud_t,aggr=None):
  return np.squeeze(P_gud_t)[:]

def get_infections(R1s,tvec,t,aggrop=None,R2s=None,vsop='1-2'):
  # TODO: upddate for foi edits
  # get fully stratified infection counts by partnership/group for plotting alluvial diagram in R
  if aggrop is None: aggrop = lambda inf: np.median(inf,axis=0)
  if isinstance(R1s,dict): R1s = [R1s]
  if isinstance(R2s,dict): R2s = [R2s]
  data = [['t','p','fs','fi','ts','ti','infections']]
  for p in range(4):
    for fs in range(2):
      for fi in range(4):
        for ts in range(2):
          for ti in range(4):
            kwds = dict(tvec=tvec,t=t,p=p,fpop=dict(s=fs,i=fi),tpop=dict(s=ts,i=ti),aggr=True)
            if R2s is None:
              inf = aggrop([whoinfectwhom(R1,**kwds) for R1 in R1s])
            else:
              inf = aggrop([vs_fun(whoinfectwhom(R1,**kwds),
                                   whoinfectwhom(R2,**kwds),vsop) for R1,R2 in zip(R1s,R2s)])
            data += [[tk,p,fs,fi,ts,ti,infk] for tk,infk in zip(t,inf)]
  return data
