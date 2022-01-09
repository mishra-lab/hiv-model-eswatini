import os
import numpy as np
from numpy import set_printoptions
from numpy import nan as NAN
from scipy.optimize import minimize
from scipy.optimize import Bounds as bounds
from pathos.multiprocessing import ProcessingPool as ppool

set_printoptions(suppress=True,linewidth=200)

_ = None

therootpath = os.path.abspath(__file__).replace(os.path.join('code','utils','__init__.py'),'')

def rootpath(*args):
  return os.path.join(therootpath,*args)

