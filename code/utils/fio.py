import os
import csv
import numpy as np
from datetime import datetime
from PyPDF2 import PdfFileMerger as pdfm

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
  return obj

def load(fname):
  # load obj from file via numpy & avoid 0-dimensional array obj for dict and other obj types
  obj = np.load(fname+'.npy',allow_pickle=True)
  if not obj.shape: # dict etc.
    return obj[()]
  return obj # array-like

def save_csv(fname,obj):
  # equivalent formats (ordered as implemented below)
  # 1. {'A':[0,1,2],'B':[3,4,5]}
  # 2. [{'A':0,'B':3},{'A':1,'B':4},{'A':2,'B':5}]
  # 3. [['A','B'],[0,3],[1,4],[2,5]]
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

def datestamp(date=None):
  # standard formatted datestamp
  if date is None: date = datetime.now()
  return date.strftime('%Y-%m-%d')

def filehash(*fnames,root=None,N=6):
  # get a unique hash based on the current state of some files
  return ''.join([os.popen('sha1sum '+os.path.join(root,fname)).read()[0:N] for fname in fnames])

def pdfmerge(ofname,fnames,rm=False):
  # concatenate pdfs
  m = pdfm()
  for fname in fnames:
    m.append(fname)
  with open(ofname,'wb') as f:
    m.write(f)
  if rm:
    for fname in fnames:
      os.remove(fname)


