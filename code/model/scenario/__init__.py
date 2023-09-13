import os
from utils import fio,rootpath,genpath,parallel
from model import system

uid = fio.datestamp()

N = dict(
  sam    = 1000,    # total initial samples
  batch  = 10,      # number of batches
  imis   = 10,      # numer of imis iterations
  fisam  = .1,      # fraction to resample per imis iter
  fpost  = .03,     # fraction to resample as posterior
)
tvec = dict(
  cal  = system.get_t(tf=2025),
  main = system.get_t(tf=2050),
  plot = system.get_t(tf=2050,dt=1),
  outs = system.get_t(tf=2050,t0=2000,dt=5),
)

def fname(ftype,phase,key,case='base',b='all'):
  subdirs = [uid,str(N['sam'])+'x'+str(N['imis'])]
  if ftype=='npy':
    path,ext = ['data','npy',*subdirs],''
  if ftype=='csv':
    path,ext = ['data','csv',*subdirs],'.csv'
  if ftype=='fig':
    path,ext = ['out','fig',*subdirs],'.pdf'
  # e.g. data/npy/2022-01-01/1000/fit_Ps_base_3.npy
  return genpath(rootpath(*path,'{}_{}_{}_{}{}'.format(phase,key,case,b,ext)))

def get_seeds(b):
  return list(range(N['bsam']*b,N['bsam']*(b+1)))

def set_N_all():
  # finalize N (dependent)
  N['bsam']  = int(N['sam']/N['batch'])  # number of samples per batch
  N['bisam'] = int(N['bsam']*N['fisam']) # number of resamples per batch-iter
  N['psam']  = int(N['sam']*N['fpost'])  # number of resamples as posterior

akwds = fio.argvkwds(scinet=False)

if akwds.pop('scinet'):
  os.environ['MPLCONFIGDIR'] = fio.tmpfile() # TODO: this may not work?
  parallel.cpus = 80

uid = '2023-09-12'
set_N_all()
