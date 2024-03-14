import numpy as np

def log(lvl,msg=''):
  # print with immediate output & hierarchy format
  if lvl < 0: print('',flush=True); return msg
  pre = ['-'*80+'\n','',' > ',''][lvl]
  end = ['\n','\n','\n',''][lvl]
  print(pre+str(msg),end=end,flush=True)

def flatten(x):
  # flatten nested lists & ensure x is iterable
  # e.g. flatten([[1,[2,[3]]],4,[[5]]]) returns [1,2,3,4,5]
  # behaviour for dicts is undefined
  f = []
  if isinstance(x,str) or not hasattr(x,'__iter__'):
    f.append(x)
  else:
    for xi in x:
      f.extend(flatten(xi))
  return f

def unique(x):
  # return unique elements in x (preserve order)
  seen = set()
  return [xi for xi in x if xi not in seen and (seen.add(xi) or True)]

def interval_qs(interval):
  # e.g. interval_qs(.95) returns (.025,.975)
  return ((1-interval)/2, 1-(1-interval)/2)

def linear_comb(eps,x1,x2):
  # e.g. linear_comb(.8,4,5) returns 4.2
  return eps * x1 + (1-eps) * x2

def clr_interp(clr1,clr2,f=.5):
  # interpolate (r,g,b) colors
  return tuple((f)*clr1[i]+(1-f)*clr2[i] for i in range(3))

def dict_split(d,keys):
  # pop keys from d if possible (modifies in-place), and return resulting dict
  return {key:d.pop(key) for key in keys if key in d}

def dict_str(d):
  # e.g. dict_str({'a':1}) returns 'a=1'
  return ', '.join([ '{}={}'.format(k,v) for k,v in d.items() ])

def dict_list_update(ds,lkwds=None,**kwds):
  # update every dict in ds with kwds and lkwds
  for i,d in enumerate(ds):
    d.update(kwds)
    if lkwds:
      d.update({k:v[i] for k,v in lkwds.items()})
  return ds

def dtfun(t):
  # compute dt from t, repeating the last dt so lengths match
  return np.diff(t,append=2*t[-1]-t[-2])

def itslice(t,tvec):
  # get boolean array for which t are in tvec
  # e.g. itslice([0,1,2,3],[0,3]) returns [True,False,False,True]
  return np.in1d(tvec,t)

def tdt(t,dt,cast=int):
  return [cast(ti) for ti in t if ti%dt==0]

def nan_to_value(x,v):
  # faster than np.nan_to_num for some reason
  x[np.isnan(x)] = v
  return x

def rk4step(Xi,ti,dt,dXfun,keys=None,**kwds):
  # Runge-Kutta 4-th order step, assuming dXfun returns a dict with at least 'dX'
  # Apply RK4 to all (default) or specified keys in the dict
  R1 = dXfun(Xi,              ti,     **kwds)
  R2 = dXfun(Xi+dt*R1['dX']/2,ti+dt/2,**kwds)
  R3 = dXfun(Xi+dt*R2['dX']/2,ti+dt/2,**kwds)
  R4 = dXfun(Xi+dt*R3['dX'],  ti+dt,  **kwds)
  if keys is None: keys = R1.keys()
  return { key: (R1[key] + 2*R2[key] + 2*R3[key] + R4[key])/6 for key in keys }

def xdi(X,di):
  # slice X using {dim:index, ...} or x.sum(axis=d) if index=None
  # never changing the number of dimensions (keepdims=True)
  # e.g. xdi(np.ones(5,4,3,2),{0:None,1:(1,2),3:0}).shape = (1,2,3,1)
  _ = None
  for d,i in di.items():
    if d==0:
      X = X.sum(axis=0,keepdims=True) if i is None \
        else X[i,_] if isinstance(i,int) \
        else X[i]
    elif d==1:
      X = X.sum(axis=1,keepdims=True) if i is None \
        else X[:,i,_] if isinstance(i,int) \
        else X[:,i]
    elif d==2:
      X = X.sum(axis=2,keepdims=True) if i is None \
        else X[:,:,i,_] if isinstance(i,int) \
        else X[:,:,i]
    elif d==3:
      X = X.sum(axis=3,keepdims=True) if i is None \
        else X[:,:,:,i,_] if isinstance(i,int) \
        else X[:,:,:,i]
    elif d==4:
      X = X.sum(axis=4,keepdims=True) if i is None \
        else X[:,:,:,:,i,_] if isinstance(i,int) \
        else X[:,:,:,:,i]
    elif d==5:
      X = X.sum(axis=5,keepdims=True) if i is None \
        else X[:,:,:,:,:,i,_] if isinstance(i,int) \
        else X[:,:,:,:,:,i]
    else:
      raise NotImplementedError('xdi @ d > 5 (d = {})'.format(d))
  return X
