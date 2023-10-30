import os
import numpy as np
from model import params
from model.scenario import fname,uid,nid
from utils import fio

def tfname(f):
  return fio.rootpath('code','params','tab',f)

def fmt(x,sf=3):
  if x > 1000: return round(x,1)
  if x == 0: return 0
  return round(x,sf-1-int(np.floor(np.log10(abs(x)))))

tex = dict(
  tab = fio.load_txt(tfname('t.tab.tex')),
  row = fio.load_txt(tfname('t.row.tex')),
)

PX = fio.load_csv(tfname('par.defs.csv'),fmt='dict') # definitions
PD = params.def_sample_distrs()                      # prior
Ps = fio.load(fname('npy','fit','Ps'))               # posterior

for X in PX:
  k = X['parameter']
  print(k)
  v = np.array([P[k] for P in Ps])
  X.update(
    imu=fmt(PD[k].mean()),
    ilo=fmt(PD[k].ppf(.025)),
    ihi=fmt(PD[k].ppf(.975)),
    omu=fmt(np.mean(v)),
    olo=fmt(np.quantile(v,.025)),
    ohi=fmt(np.quantile(v,.975)),
    distr=PD[k].dist.name.capitalize())
rows = ''.join(tex['row'].format(**X) for X in PX)

# NOTE: relies on docs/app/config.tex
fio.save_txt(fio.rootpath('out','tex',uid,nid,'tab.par.tex'),
  tex['tab'].replace('{{ rows }}',rows))
