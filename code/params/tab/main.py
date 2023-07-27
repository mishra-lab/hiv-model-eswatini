import numpy as np
from model import params
from model.scenario import fname
from utils import fio

def csvname(key):
  return fio.rootpath('code','params','tab','pars_{}.csv'.format(key))

def fmt(x,sf=3):
  if x > 1000: return round(x,1)
  return round(x,sf-1-int(np.floor(np.log10(abs(x)))))

PD = params.def_sample_distrs()               # prior
Ps = fio.load(fname('npy','fit','Ps'))        # posterior
PX = fio.load_csv(csvname('def'),fmt='rows')  # definitions

PX[0] += ['imu','ilo','ihi', # prior:     mean, 95.lo, 95.hi
          'omu','olo','ohi'] # posterior: mean, 95.lo, 95.hi
for row in PX[1:]:
  key = row[0]
  dk,xk = PD[key],[P[key] for P in Ps]
  row += [fmt(x) for x in [
    dk.mean(),*dk.interval(.95),
    np.mean(xk),*np.quantile(xk,(.025,.975))]]

fio.save_csv(csvname('full'),PX)
