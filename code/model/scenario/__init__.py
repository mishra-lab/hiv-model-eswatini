import os
from utils import log,fio,rootpath,genpath,parallel
from model import system

uid = fio.datestamp()

N = dict(
  batch  =  100,    # [ 100] number of batches
  hsam   = 1000,    # [1000] initial hypercube samples per batch
  isam   =  100,    # [ 100] number to resample per imis iter
  imis   =  100,    # [ 100] number of imis iterations
  post   = 1000,    # [1000] number of posterior samples
) # h1000i100b100
# N = dict(batch=5,hsam=100,isam=10,imis=15,post=100) # DEBUG: h100i15b5
tvec = dict(
  cal  = system.get_t(tf=2025),
  main = system.get_t(tf=2050),
  plot = system.get_t(tf=2050,dt=1),
  outs = system.get_t(tf=2050,t0=2000,dt=5),
)

def fname(ftype,phase,key,case='base',b='all'):
  nid = 'h{hsam}i{imis}b{batch}'.format(**N)
  subdirs = [uid,nid]
  if ftype=='npy':
    path,ext = ['data','npy',*subdirs],'.npy'
  if ftype=='csv':
    path,ext = ['data','csv',*subdirs],'.csv'
  if ftype=='fig':
    path,ext = ['out','fig',*subdirs],'.pdf'
  # e.g. data/npy/2023-01-01/h1000i25b50/fit_Ps_base_3.npy
  return genpath(rootpath(*path,'{}_{}_{}_{}{}'.format(phase,key,case,b,ext)))

def get_seeds(b):
  return list(range(N['hsam']*b,N['hsam']*(b+1)))

uid = '2024-03-12'
nid = 'h{hsam}i{imis}b{batch}'.format(**N)
log(0,'{}/{}'.format(uid,nid))
