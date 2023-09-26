import os
from utils import fio,rootpath,genpath,parallel
from model import system

uid = fio.datestamp()

N = dict(
  batch  = 1000,    # [1000] number of batches
  bsam   = 1000,    # [1000] initial samples per batch
  isam   = 100,     # [ 100] number to resample per imis iter
  imis   = 25,      # [  25] numer of imis iterations
)
tvec = dict(
  cal  = system.get_t(tf=2025),
  main = system.get_t(tf=2050),
  plot = system.get_t(tf=2050,dt=1),
  outs = system.get_t(tf=2050,t0=2000,dt=5),
)

def fname(ftype,phase,key,case='base',b='all'):
  subdirs = [uid,str(N['bsam'])+'x'+str(N['imis'])]
  if ftype=='npy':
    path,ext = ['data','npy',*subdirs],''
  if ftype=='csv':
    path,ext = ['data','csv',*subdirs],'.csv'
  if ftype=='fig':
    path,ext = ['out','fig',*subdirs],'.pdf'
  # e.g. data/npy/2023-01-01/1000x25/fit_Ps_base_3.npy
  return genpath(rootpath(*path,'{}_{}_{}_{}{}'.format(phase,key,case,b,ext)))

def get_seeds(b):
  return list(range(N['bsam']*b,N['bsam']*(b+1)))

akwds = fio.argvkwds(scinet=False)

if akwds.pop('scinet'):
  os.environ['MPLCONFIGDIR'] = fio.tmpfile() # TODO: this may not work?
  parallel.cpus = 80

uid = '2023-09-20'
