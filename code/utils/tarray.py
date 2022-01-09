from utils import deco
import numpy as np

# TODO: slicing?

class tarray:
  def __init__(self,ti,xi):
    ti,xi = np.array(ti),np.array(xi)
    self.shape = xi.shape[0:-1]
    self.params = self.fit(ti,xi)

  def __call__(self,t):
    tsize = np.array(t).size
    if tsize == 1:
      return np.reshape([eval_spline(t,**p) for p in self.params],self.shape)
    else:
      return np.reshape([eval_spline(ti,**p) for p in self.params for ti in t],(*self.shape,tsize))

  def fit(self,ti,xi):
    return [fit_spline(ti,xi[(*i,slice(None))]) for i in np.ndindex(self.shape)]

  def reshape(self,shape):
    size = np.prod(self.shape)
    if size != np.prod(shape):
      raise ValueError('cannot reshape tarray of size {} into shape {}'.format(size,shape))
    self.shape = shape
    return self

@deco.nowarn
def fit_spline(ti,xi):
  # monotonic
  i = np.isfinite(xi)
  ti,xi = ti[i],xi[i]
  n = len(ti)
  dti = np.diff(ti)
  dxi = np.diff(xi)
  di = dxi/dti
  qi = dti[0:n-2]+dti[1:n-1]
  c1 = 3*qi / ((qi+dti[1:n-1])/di[0:n-2] + (qi+dti[0:n-2])/di[1:n-1])
  c1 = np.array([di[0],*c1,di[-1]])
  qi = c1[0:n-1] + c1[1:n] - 2*di[0:n-1]
  c2 = (di-c1[0:n-1]-qi)/dti
  c2 = np.array([*c2,c2[-1]])
  c3 = qi/dti/dti
  c3 = np.array([*c3,c3[-1]])
  return dict(ti=ti,xi=xi,c1=c1,c2=c2,c3=c3)

def eval_spline(t,ti,xi,c1,c2,c3):
  i = ti.searchsorted(t)-1
  i = np.clip(i,0,len(ti)-2)
  dt = t - ti[i]
  return xi[i] + c1[i]*dt + c2[i]*dt**2 + c3[i]*dt**3
