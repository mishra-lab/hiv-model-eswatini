import os
import csv
import numpy as np
from datetime import datetime
from PyPDF2 import PdfFileMerger as pdfm
from utils.log import log

therootpath = os.path.abspath(__file__).replace(os.path.join('code','utils','fio.py'),'')

def rootpath(*args):
  return os.path.join(therootpath,*args)

def genpath(fname):
  path = os.path.dirname(fname)
  if path: os.makedirs(path,exist_ok=True)
  return fname

def save(fname,obj):
  # save obj to file via numpy, ensuring the path exists
  log(2,'fio.save: '+fname)
  np.save(genpath(fname),obj)
  return obj

def load(fname):
  # load obj from file via numpy & avoid 0-dimensional array obj for dict and other obj types
  log(2,'fio.load: '+fname)
  obj = np.load(fname+'.npy',allow_pickle=True)
  if not obj.shape: # dict etc.
    return obj[()]
  return obj # array-like

def save_csv(fname,obj):
  # equivalent formats (ordered as implemented below)
  # 1. {'A':[0,1,2],'B':[3,4,5]}
  # 2. [{'A':0,'B':3},{'A':1,'B':4},{'A':2,'B':5}]
  # 3. [['A','B'],[0,3],[1,4],[2,5]]
  log(2,'fio.save_csv: '+fname)
  with open(genpath(fname),'w') as f:
    if isinstance(obj,dict):
      w = csv.writer(f)
      w.writerow(obj.keys())
      w.writerows(zip(*obj.values()))
    elif isinstance(obj[0],dict):
      w = csv.DictWriter(f,obj[0].keys())
      w.writeheader()
      w.writerows(obj)
    else:
      w = csv.writer(f)
      w.writerows(obj)

def load_csv(fname,fmt='rows',cast=float):
  # fmt options are same as save_csv (default=2)
  # 1. fmt='cols': {'A':[0,1,2],'B':[3,4,5]}
  # 2. fmt='dict': [{'A':0,'B':3},{'A':1,'B':4},{'A':2,'B':5}]
  # 3. fmt='rows': [['A','B'],[0,3],[1,4],[2,5]]
  log(2,'fio.load_csv: '+fname)
  if cast is None: case = lambda x: x
  with open(fname,'r') as f:
    if fmt=='cols':
      r = csv.reader(f)
      keys = next(r)
      return {key:[cast(x) for x in col] for key,col in zip(keys,list(zip(*r)))}
    if fmt=='dict':
      r = csv.DictReader(f)
      return [{key:cast(x) for key,x in row.items()} for row in r]
    if fmt=='rows':
      r = csv.reader(f)
      keys = next(r)
      return [keys]+[[cast(x) for x in row] for row in r]

def datestamp(date=None):
  # standard formatted datestamp
  if date is None: date = datetime.now()
  return date.strftime('%Y-%m-%d')

def filehash(*fnames,root=None,n=7):
  # get a unique hash based on the current state of some files
  root = '.' if root is None else root
  return ''.join([os.popen('sha1sum '+os.path.join(root,fname)).read()[0:n] for fname in fnames])

def randhash(n=7):
  # get a random hash
  chars = '0123456789abcdefghijklmnopqrstuvwxyz'
  return ''.join(np.random.choice(list(chars),n))

def tmpfile(fname=None,root=None,n=7):
  # generate a dir based on randhash & return path
  fname = '' if fname is None else fname
  root  = '.tmp' if root is None else root
  return genpath(os.path.join(root,randhash(n),fname))

def pdfmerge(ofname,fnames,rm=False):
  # concatenate pdfs
  log(2,'fio.pdfmerge: '+ofname)
  m = pdfm()
  for fname in fnames:
    m.append(fname)
  with open(ofname,'wb') as f:
    m.write(f)
  if rm:
    for fname in fnames:
      os.remove(fname)


