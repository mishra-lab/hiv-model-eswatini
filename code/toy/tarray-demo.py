import numpy as np
from numpy import nan
from utils import tarray as ta
import matplotlib.pyplot as plt

t = np.arange(2000,2030,.1)
to = [1999,2000,2005,2015,2030,2031]
Xo = np.array([
     [ .00, .00, .35, .55, .90, .90],
     [ .00, .00, nan, .30, .40, .40],
  ])
X = ta.tarray(to,Xo)

clr = 'rb'
for i in range(Xo.shape[0]):
  for ti,xij in zip(to,Xo[i]):
    plt.plot(ti,xij,clr[i]+'.')
  plt.plot(t,X(t)[i],clr[i]+'-',alpha=.5)
plt.show()
