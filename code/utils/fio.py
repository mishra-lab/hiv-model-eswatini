import os
import numpy as np

therootpath = os.path.abspath(__file__).replace(os.path.join('code','utils','fio.py'),'')

def rootpath(*args):
  return os.path.join(therootpath,*args)

def genpath(fname):
  path = os.path.dirname(fname)
  if path: os.makedirs(path,exist_ok=True)
  return fname

def save(fname,obj):
  # save obj to file via numpy, ensuring the path exists
  np.save(genpath(fname),obj)

def load(fname):
  # load obj from file via numpy & avoid 0-dimensional array obj for dict and other obj types
  obj = np.load(fname+'.npy',allow_pickle=True)
  if not obj.shape: # dict etc.
    return obj[()]
  return obj # array-like

def filehash(*fnames,root=None,N=6):
  # get a unique hash based on the current state of some files
  return ''.join([os.popen('sha1sum '+os.path.join(root,fname)).read()[0:N] for fname in fnames])
