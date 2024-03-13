import os
import sys
import csv
import numpy as np
from datetime import datetime
from PyPDF2 import PdfMerger as pdfm
from utils.log import log

# therootpath: full path to parent of /code/utils/fio.py
# so we can reliably create full paths regardless of the project root folder location
therootpath = os.path.abspath(__file__).replace(os.path.join('code','utils','fio.py'),'')

def rootpath(*args):
  # e.g. rootpath('a','b') returns therootpath/a/b
  return os.path.join(therootpath,*args)

def genpath(fname):
  # create path to fname if needed
  path = os.path.dirname(fname)
  if path: os.makedirs(path,exist_ok=True)
  return fname

def save_npy(fname,obj):
  # save obj to file via numpy, ensuring the path exists
  log(2,'fio.save_npy: '+fname)
  np.save(genpath(fname),obj)
  return obj

def load_npy(fname):
  # load obj from file via numpy & avoid 0-dimensional array obj for dict and other obj types
  log(2,'fio.load_npy: '+fname)
  obj = np.load(fname,allow_pickle=True)
  if not obj.shape: # dict etc.
    return obj[()]
  return obj # array-like

def save_csv(fname,obj):
  # equivalent formats (ordered as implemented below)
  # 1. cols: {'A':[0,1,2],'B':[3,4,5]}
  # 2. dict: [{'A':0,'B':3},{'A':1,'B':4},{'A':2,'B':5}]
  # 3. rows: [['A','B'],[0,3],[1,4],[2,5]]
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

def load_csv(fname,fmt='rows',cast=float,castbak=str,**kwds):
  # fmt options are same as save_csv
  # 1. 'cols': {'A':[0,1,2],'B':[3,4,5]}
  # 2. 'dict': [{'A':0,'B':3},{'A':1,'B':4},{'A':2,'B':5}]
  # 3. 'rows': [['A','B'],[0,3],[1,4],[2,5]]
  # we try to cast variables using 'cast', and default to castbak if it fails
  log(2,'fio.load_csv: '+fname)
  if cast is None: cast = lambda x: x
  def castfun(x):
    try: return cast(x)
    except: return castbak(x)
  kwds = dict(**kwds,skipinitialspace=True)
  with open(fname,'r') as f:
    if fmt=='cols':
      r = csv.reader(f,**kwds)
      keys = next(r)
      return {key:[castfun(x) for x in col] for key,col in zip(keys,list(zip(*r)))}
    if fmt=='dict':
      r = csv.DictReader(f,**kwds)
      return [{key:castfun(x) for key,x in row.items()} for row in r]
    if fmt=='rows':
      r = csv.reader(f,**kwds)
      keys = next(r)
      return [keys]+[[castfun(x) for x in row] for row in r]

def load_txt(fname):
  log(2,'fio.load_txt: '+fname)
  with open(fname,'r') as f:
    return(f.read())

def save_txt(fname,obj):
  log(2,'fio.save_txt: '+fname)
  with open(fname,'w') as f:
    return(f.write(str(obj)))

def datestamp(date=None):
  # standard formatted datestamp
  if date is None: date = datetime.now()
  return date.strftime('%Y-%m-%d')

def filehash(*fnames,root=None,n=7):
  # get a unique n-length hash based on the current state of some file(s)
  root = '.' if root is None else root
  return ''.join([os.popen('sha1sum '+os.path.join(root,fname)).read()[0:n] for fname in fnames])

def randhash(n=7):
  # get a random n-length hash
  chars = '0123456789abcdefghijklmnopqrstuvwxyz'
  return ''.join(np.random.choice(list(chars),n))

def tmpfile(fname=None,root=None,n=7):
  # generate a dir based on randhash, create, & return path
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

def argvkwds(i=1,**kwds): # delete
  # parse argv to dict, e.g. ['a','b=1'] -> dict(a=True,b=1)
  for arg in sys.argv[i:]:
    k,e,v = arg.partition('=')
    if e:
      try:
        kwds[k] = int(v)
      except:
        kwds[k] = v
    else:
      kwds[k] = True
  return kwds
