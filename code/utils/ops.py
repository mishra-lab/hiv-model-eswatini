import numpy as np
from utils import os,rootpath
from scipy.optimize import minimize
from scipy.optimize import Bounds as bounds
from pathos.multiprocessing import ProcessingPool as ppool

# TODO: should this just go in __init__.py ?

_glob = {} # see globs

def flatten(x):
  # flatten & ensure x is iterable
  # e.g. [[1,[2,[3]]],4,[[5]]] -> [1,2,3,4,5]
  f = []
  if isinstance(x,str) or not hasattr(x,'__iter__'):
    f.append(x)
  else:
    for xi in x:
      f.extend(flatten(xi))
  return f

def around(x,decimals):
  # TODO: remove?
  try:
    return np.around(x,decimals)
  except:
    return x

def keystartswith(d,s):
  # convenience: get all keys in "d" starting with "s"
  return [key for key in d.keys() if key.startswith(s)]

def intervalqs(interval):
  return ((1-interval)/2, 1-(1-interval)/2)

def lincomb(eps,x1,x2):
  # convenience: see code
  return eps * x1 + (1-eps) * x2

def cinterp(c1,c2,f=.5):
  return tuple((f)*c1[i]+(1-f)*c2[i] for i in range(3))

def dictsplit(d,keys):
  # pop keys from d if possible, and return resulting dict
  return {key:d.pop(key) for key in keys if key in d}

def dictstr(d,decimals=4):
  # print a dictionary "pretty"
  return ', '.join([
    '{}={}'.format(k,around(v,decimals)) for k,v in d.items()
  ])

def squarish(n):
  # how many rows,cols sould we use to arrange subplots, given n
  # usually "squarish", except for some special cases (fixed)
  col = int(np.ceil(np.sqrt(n)))
  row = int(np.ceil(n/col))
  fixed = {3:(1,3),4:(1,4),7:(2,4),8:(2,4),32:(4,8)}
  return fixed.get(n,(row,col))

def npsave(fname,obj):
  # save obj to file via numpy, ensuring the path exists
  rootfname = rootpath('data','.npy',fname)
  os.makedirs(os.path.dirname(rootfname),exist_ok=True)
  np.save(rootfname,obj)

def npload(fname):
  # load obj from file via numpy & avoid 0-dimensional array obj for dict and other obj types
  obj = np.load(rootpath('data','.npy',fname+'.npy'),allow_pickle=True)
  if not obj.shape: # dict etc.
    return obj[()]
  return obj # array-like

def filehash(*fnames,root=None,N=6):
  # get a unique hash based on the current state of some files
  return ''.join([os.popen('sha1sum '+os.path.join(root,fname)).read()[0:N] for fname in fnames])

def globs(get=None,clear=None,iter=False,default=None,**kwds):
  # tool for getting, clearing, and setting global variables (one dictionary)
  # get:     key to get
  # clear:   key to clear
  # **kwds:  "key=value" variable assignments
  # iter:    for kwds: is value iterable & should we append to it, vs overwrite?
  # default: for get & clear: default get value if ket not found
  global _glob
  if get:
    return _glob.get(get,default)
  if clear:
    return _glob.pop(clear,default)
  for key,item in kwds.items():
    if iter:
      if key in _glob:
        _glob[key].append(item)
      else:
        _glob[key] = [item]
    else:
      _glob[key] = item
  return None