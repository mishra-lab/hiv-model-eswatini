import numpy as np

def monotonic(x,axis=0,sign=1):
  x = np.moveaxis(x,axis,0)
  x2 = x.reshape((x.shape[0],-1))
  v2 = np.isfinite(x2)
  for i in range(x2.shape[1]):
    if any(np.sign(np.diff(x2[v2[:,i],i])) != sign):
      return False
  return True

def dimarray(dims):
  return np.meshgrid(*dims.values(),indexing='ij')

def melt(X,dims):
  return([list(dx) for dx in
    zip(*[Mi.flatten() for Mi in [*dimarray(dims),X]])
  ])

def meltstr(X,dims,fmt='{:9.3f}'):
  D = [*dimarray(dims)]
  A = D[0]
  for Di in D[1:]:
    A = np.char.add(A,',')
    A = np.char.add(A,Di)
  f = np.vectorize(lambda v: fmt.format(v))
  A = np.char.add(A,': ')
  A = np.char.add(A,f(X))
  return(A)

def dtfun(t):
  return np.diff(t,append=2*t[-1]-t[-2])