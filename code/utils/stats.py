import scipy.stats as ss
from scipy.stats.qmc import LatinHypercube
from numpy.random import permutation
import numpy as np

def beta_binom(p,n):
  # beta approximation of binomial
  return ss.beta(a=p*n, b=(1-p)*n)

def gamma_p(p,v):
  # gamma distribution using 
  return ss.gamma(a=p*p/v,scale=v/p)

def uniform(l,h):
  return ss.uniform(loc=l,scale=h-l)

def ratio_binom(p1,n1,p2,n2):
  # Katz1978
  var = ((1/p1)-1)/n1+((1/p2)-1)/n2
  return ss.lognorm(scale=p1/p2,s=np.sqrt(var))

def lhs(d,n,seed=None):
  # n: number of samples to obtain
  return LatinHypercube(d,seed=seed).random(n)
