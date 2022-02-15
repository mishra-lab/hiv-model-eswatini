import sys
from utils import fio,parallel,linear_comb
from model import system,params,target,fit,scenario,out
import numpy as np

scenario.uid = '2022-02-13'
scenario.N['size'] = 70 # DEBUG

scenario.main_fit(sample=True,refit=True)
scenario.main_random(sample=True,resample=True)
