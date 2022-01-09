import warnings
import numpy as np

def nowarn(fun):
  # don't print warnings while executing fun
  def decorator(*args,**kwds):
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      result = fun(*args,**kwds)
      return result
  return decorator

def rmap(args=['X']):
  # pass values from "R" as arguments to fun
  def wrapper(fun):
    def decorator(R,**kwds):
      kwds.update({k:R[k] for k in args})
      return fun(**kwds)
    return decorator
  return wrapper

def tslice(targs=[]):
  # slice args with t before computing fun, assuming first dimension is t dimension
  # targs: list of indices for args that need slicing (fixed for fun)
  # t:     t values that fun should compute for
  # tvec:  complete list of t values (for look-up)
  def wrapper(fun):
    def decorator(t=None,tvec=None,**kwds):
      if t is not None:
        it = np.in1d(tvec,t) # boolean
        for k in kwds.keys():
          if k in targs:
            kwds[k] = kwds[k][it] # slicing (keeps singleton)
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
