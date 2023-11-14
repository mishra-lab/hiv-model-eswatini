import numpy as np
from itertools import product as iprod
from utils import _,deco,xdi,dtfun,itslice
from model import system,foi,slicers
# TODO: integrate slicers fully?

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

@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def NX(X,s=None,i=None,aggr=True):
  X = X.sum(axis=(3,4))
  X = xdi(X,{1:s,2:i})
  return X.sum(axis=(1,2)) if aggr else X

@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def Psi(X,s=None,i=None,aggr=True):
  X  = X.sum(axis=(3,4))
  XS = X.sum(axis=(1,2),keepdims=True)
  X  = xdi(X,{1:s,2:i})
  return aggratio(X,XS,aggr)

@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def prevalence(X,s=None,i=None,aggr=True):
  X  = xdi(X,{1:s,2:i})
  XS = X.sum(axis=(3,4)) # <- denominator; numerator -> (1 - susceptible)
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4))
  return aggratio(Xhiv,XS,aggr)

@deco.rmap(Rk=['X','inc'],Pk=['foi_mode'])
@deco.tslice(tk=['X','inc'])
def incidence(X,inc,foi_mode,s=None,i=None,aggr=True):
  inf = foi.aggr_inc(inc,foi_mode,axis=(1,4,5),Xsus=X[:,:,:,0,0])
  inf_si = xdi(inf,{1:s,2:i})
  sus_si = xdi(X[:,:,:,0,0],{1:s,2:i})
  return aggratio(inf_si,sus_si,aggr)

@deco.rmap(Rk=['X','inc'],Pk=['foi_mode'])
def cuminfect(X,inc,foi_mode,tvec,s=None,i=None,aggr=True,t0=None):
  dt = dtfun(tvec)
  inf = foi.aggr_inc(inc,foi_mode,axis=(1,4,5),Xsus=X[:,:,:,0,0])
  inf_si = xdi(inf,{1:s,2:i})
  inf_dt = inf_si.sum(axis=(1,2)) * dt if aggr else inf_si * dt[:,_,_]
  if t0:
    inf_dt[tvec < t0] = 0
  return np.cumsum(inf_dt,axis=0)

@deco.nanzero
@deco.rmap(Rk=['Xk'],Pk=['K_psi'])
@deco.tslice(tk=['Xk'])
def tdsc(Xk,K_psi,p=None,s=None,i=None,aggr=True,sus=False):
  # NOTE: only works for foi_mode='base' & system.run(...,Xk=True)
  # Xk.shape: (t:*, s:2, i:4, k:5, h:6, c:5) -> (t:*, k:5, s:2, i:4)
  Xk = np.moveaxis(Xk if sus else Xk[:,:,:,:,1:,:],3,1).sum(axis=(4,5)) # only inf or sus too
  XKsc = xdi(Xk[:,1:,:,:],{1:p,2:s,3:i})
  XK   = xdi(xdi(Xk,{1:_}) * np.squeeze(K_psi)[_,:,:,:],{1:p,2:s,3:i})
  return aggratio(XKsc,XK,aggr,axis=(1,2,3))

def can_tdsc(R):
  # check if we can compute tdsc from this R
  return ('Xk' in R) and (R['P']['foi_mode'] == 'base')

@deco.nanzero
@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def Ph(X,h,s=None,i=None,aggr=True):
  X    = xdi(X,{1:s,2:i})
  Xh   = xdi(X,{3:h}).sum(axis=(3,4))
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4))
  return aggratio(Xh,Xhiv,aggr)

@deco.nanzero
@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def diagnosed(X,s=None,i=None,aggr=True):
  X = xdi(X,{1:s,2:i})
  Xhiv = X[:,:,:,1:,:].sum(axis=(3,4)) # PLHIV
  Xdia = X[:,:,:,1:,1:].sum(axis=(3,4)) # diagnosed
  return aggratio(Xdia,Xhiv,aggr)

@deco.nanzero
@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def treated(X,s=None,i=None,aggr=True,cond=False):
  X = xdi(X,{1:s,2:i})
  if cond:
    Xref = X[:,:,:,1:,1:].sum(axis=(3,4)) # diagnosed
  else:
    Xref = X[:,:,:,1:,:].sum(axis=(3,4)) # PLHIV
  Xtre = X[:,:,:,1:,3:].sum(axis=(3,4))
  return aggratio(Xtre,Xref,aggr)

@deco.nanzero
@deco.rmap(Rk=['X'])
@deco.tslice(tk=['X'])
def vls(X,s=None,i=None,aggr=True,cond=False):
  X = xdi(X,{1:s,2:i})
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
@deco.rmap(Rk=['X','dx_sit'])
@deco.tslice(tk=['X','dx_sit'])
def dx_rate(X,dx_sit,s=None,i=None,aggr=True):
  X = X[:,:,:,1:,0].sum(axis=(3)) # undiagnosed only
  Xdx = xdi(X*np.squeeze(dx_sit),{1:s,2:i})
  XS  = xdi(X,{1:s,2:i})
  return aggratio(Xdx,XS,aggr)

@deco.nanzero
@deco.rmap(Rk=['X','tx_sit','Rtx_ht'])
@deco.tslice(tk=['X','tx_sit','Rtx_ht'])
def tx_rate(X,tx_sit,Rtx_ht,s=None,i=None,aggr=True):
  X = X[:,:,:,1:,1] # diagnosed only
  Xtx = xdi(X*np.squeeze(tx_sit*Rtx_ht),{1:s,2:i}).sum(axis=3)
  XS  = xdi(X,{1:s,2:i}).sum(axis=3)
  return aggratio(Xtx,XS,aggr)

@deco.rmap(Rk=['PF_condom_t'])
@deco.tslice(tk=['PF_condom_t'])
def condom(PF_condom_t,p,aggr=None):
  return np.squeeze(PF_condom_t)[:,p]

@deco.rmap(Rk=['PF_circum_t'])
@deco.tslice(tk=['PF_circum_t'])
def circum(PF_circum_t,aggr=None):
  return np.squeeze(PF_circum_t)[:]

@deco.rmap(Rk=['X','inc'],Pk=['foi_mode'])
@deco.tslice(tk=['X','inc'])
def infections(X,inc,foi_mode,p,fs,fi,ts,ti):
  # NOTE: p,fs,fi,ts,ti must be single values!
  return foi.aggr_inc(inc[:,p,ts,ti,fs,fi],foi_mode,axis=(),Xsus=X[:,ts,ti,0,0])

def wiw(R1s,tvec,t,R2s=None,vsop='1-2'):
  if isinstance(R1s,dict): R1s = [R1s]
  if isinstance(R2s,dict): R2s = [R2s]
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

def expo(R1s,tvec,t,onames,snames,R2s=None,vsop='raw',ecols=None,mode='q',**kwds):
  # E: dict of cols (lists); first cols = strata, later cols = output quantiles/values per id
  if mode == 'q':
    aggrop = lambda os: np.nanquantile(os,qs,axis=0)
    cols = ['q'+str(q) for q in qs]
  if mode == 'id':
    aggrop = lambda os: np.array(os)
    cols = ['i'+str(R['P']['id']) for R in R1s]
  sg,og,tg = [g.flatten().tolist() for g in np.meshgrid(snames,onames,t)]
  if ecols is None: ecols = {}
  ecols.update(op=vsop)
  E = dict(out=og,pop=sg,t=tg,**{k:[v]*len(tg) for k,v in ecols.items()},**{k:[] for k in cols})
  for oname in onames:
    fun = by_name(oname)
    for sname in snames:
      if oname == 'cuminfect': # special case as we cannot use @deco.tslice
        sfun = lambda R: fun(R,**slicers[sname].pop,tvec=tvec,**kwds)[itslice(t,tvec)]
      else:
        sfun = lambda R: fun(R,**slicers[sname].pop,tvec=tvec,t=t,**kwds)
      if R2s is None:
        osc = aggrop([sfun(R) for R in R1s])
      else:
        osc = aggrop([vs_fun(sfun(R1),sfun(R2),vsop) for R1,R2 in zip(R1s,R2s)])
      for i,col in enumerate(cols): # TODO: is this slow? if so, vectorize with np?
        E[col] += osc[i,:].tolist()
  return E
