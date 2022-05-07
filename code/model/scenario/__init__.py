from utils import fio,rootpath,genpath
from model import system

uid = fio.datestamp()

N = dict(
  sam    = 100000, # total samples
  topcal = .1,     # top fraction to save from each cal batch
  topfit = .01,    # top fraction to save & use from all batches
  batch  = 10,     # how many batches
)
tvec = dict(
  cal  = system.get_t(tf=2021),
  main = system.get_t(tf=2050),
  plot = system.get_t(tf=2050,dt=1),
  outs = system.get_t(tf=2050,t0=2000,dt=5),
)

def fname(ftype,phase,key,case='base',b='all'):
  subdirs = [uid,str(N['sam'])]
  if ftype=='npy':
    path,ext = ['data','npy',*subdirs],''
  if ftype=='csv':
    path,ext = ['data','csv',*subdirs],'.csv'
  if ftype=='fig':
    path,ext = ['out','fig',*subdirs],'.pdf'
  # e.g. data/npy/2022-01-01/1000/fit_Ps_base_3.npy
  return genpath(rootpath(*path,'{}_{}_{}_{}{}'.format(phase,key,case,b,ext)))

def batch_select(objs):
  nb = int(len(objs) / N['batch'])
  return objs[slice(nb*N['b'],nb*(N['b']+1))]


