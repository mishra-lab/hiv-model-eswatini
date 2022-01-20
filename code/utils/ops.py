import numpy as np

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

def unique(x):
  # return unique elements in x (preserve order)
  seen = set()
  return [xi for xi in x if xi not in seen and (seen.add(xi) or True)]

def interval_qs(interval):
  return ((1-interval)/2, 1-(1-interval)/2)

def linear_comb(eps,x1,x2):
  # convenience: see code
  return eps * x1 + (1-eps) * x2

def clr_interp(clr1,clr2,f=.5):
  # interpolate colors
  return tuple((f)*clr1[i]+(1-f)*clr2[i] for i in range(3))

def dict_split(d,keys):
  # pop keys from d if possible, and return resulting dict
  return {key:d.pop(key) for key in keys if key in d}

def dict_str(d):
  # print a dictionary "pretty"
  return ', '.join([ '{}={}'.format(k,v) for k,v in d.items() ])

def squarish(n):
  # how many rows,cols sould we use to arrange subplots, given n
  # usually "squarish", except for some special cases (fixed)
  col = int(np.ceil(np.sqrt(n)))
  row = int(np.ceil(n/col))
  fixed = {3:(1,3),4:(1,4),7:(2,4),8:(2,4),32:(4,8)}
  return fixed.get(n,(row,col))

def dtfun(t):
  return np.diff(t,append=2*t[-1]-t[-2])

def itslice(t,tvec):
  return np.in1d(tvec,t)

def nan_to_value(x,v):
  # faster than np.nan_to_num as we don't deal with infs?
  x[np.isnan(x)] = v
  return x
