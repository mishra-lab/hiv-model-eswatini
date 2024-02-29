import warnings
import numpy as np
from utils import itslice

def nowarn(fun):
  # don't print warnings while executing fun
  def decorator(*args,**kwds):
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      result = fun(*args,**kwds)
      return result
  return decorator

def rmap(Rk=[],Pk=[]):
  # pass values from 'R' and/or R['P'] as kwds to fun
  # e.g. if we define @rmap(Rk=['X']) def fun(X,a):
  # we can use it like fun(R,a) but it acts like fun(a,X=R['X'])
  # this allows us to call all out functions like out(R,...)
  def wrapper(fun):
    def decorator(R,**kwds):
      kwds.update({k:R[k] for k in Rk})
      kwds.update({k:R['P'][k] for k in Pk})
      return fun(**kwds)
    return decorator
  return wrapper

def tslice(tk=[]):
  # slice args with t before computing fun, assuming first dimension is t dimension
  # tk:   names of kwds that need slicing (always same kwds for a given fun)
  # t:    t values that fun should compute for
  # tvec: complete list of t values (for look-up)
  def wrapper(fun):
    def decorator(t=None,tvec=None,**kwds):
      if t is not None:
        it = itslice(t,tvec) # boolean
        for k in kwds.keys():
          if k in tk:
            kwds[k] = kwds[k][it] # slice (keeps singleton)
      # else: don't slice -> kwds unchanged
      return fun(**kwds)
    return decorator
  return wrapper

def nanzero(fun):
  # replace nans with 0 in the return of fun
  def decorator(*args,**kwds):
    result = nowarn(fun)(*args,**kwds)
    result[np.isnan(result)] = 0
    return result
  return decorator

def runtime(how='print'):
  # estimate the runtime of fun
  def wrapper(fun):
    def decorator(*args,**kwds):
      import time
      timer = time.time()
      result = fun(*args,**kwds)
      dt = time.time() - timer
      if how == 'return':
        return(dt)
      if how == 'print':
        print(dt)
      if how == 'dict':
        result['runtime'] = dt
      return result
    return decorator
  return wrapper
