import time
import scipy.stats as ss
from scipy.stats.qmc import LatinHypercube
from numpy.random import permutation
import numpy as np

def betabin(p,n):
  # beta approximation of binomial
  return ss.beta(a=p*n, b=(1-p)*n)

def gamma(m,sd):
  # gamma distribution using mean & sd
  return ss.gamma(a=m*m/sd/sd,scale=sd*sd/m)

def lnorm(m,sd,loc=0):
  # lognormal distribution
  q = 1+sd*sd/m/m
  return ss.lognorm(s=np.sqrt(np.log(q)),scale=m/np.sqrt(q),loc=loc)

def skewnorm(m,sd,a=0):
  d = a/np.sqrt(1+a*a)
  w = sd/np.sqrt(1-2*d*d/np.pi)
  return ss.skewnorm(a=a,loc=m-w*d*np.sqrt(2/np.pi),scale=w)

def invgauss(m,sd,z=0):
  w = m-z
  s = w**3/sd**2
  return ss.invgauss(mu=w/s,loc=z,scale=s)

def uniform(l,h):
  return ss.uniform(loc=l,scale=h-l)

def mvn(m,cov):
  # TODO: singular safe?
  return ss.multivariate_normal(mean=m,cov=cov,allow_singular=True)

def ratio_binom(p1,n1,p2,n2):
  # Katz1978
  var = ((1/p1)-1)/n1+((1/p2)-1)/n2
  return ss.lognorm(scale=p1/p2,s=np.sqrt(var))

def lhs(d,n,seed=None):
  # n: number of samples to obtain
  return LatinHypercube(d,seed=seed).random(n)

def random_seed(imax=4294967295):
  return int(imax*(time.time()%1))
