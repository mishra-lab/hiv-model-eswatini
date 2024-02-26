import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from model.target import filter_targets
from utils import _,log,genpath,flatten,dict_split,interval_qs,clr_interp,squarish,itslice,tdt
from model import out,strats

# ------------------------------------------------------------------------------
# helper functions
# TODO: some of these should go in utils (?)

def cmap(N,trim=False,option='inferno'):
  if trim: N += 2
  cmfun = cm.get_cmap(option,N)
  cmlist = [cmfun(i) for i in range(N)]
  return cmlist[1:-1] if trim else cmlist

def subplots(row,col):
  fh,ah = plt.subplots(row,col)
  ah = np.reshape(ah,(row,col))
  return fh,ah

def save(fname,*args,close=True,**kwds):
  log(2,'plot.save: '+fname)
  plt.savefig(genpath(fname),*args,**kwds)
  if close: plt.close()
  return fname

def labels(x='Year',y=None,title=None,fs=11):
  if x     is not None: plt.xlabel(x,fontsize=fs)
  if y     is not None: plt.ylabel(y,fontsize=fs)
  if title is not None: plt.title(title,fontsize=fs)

def lims(x=None,y=None):
  if x is not None: plt.xlim(x)
  if y is not None: plt.ylim(y)

# ------------------------------------------------------------------------------
# main plot functions

def line(t,x,taxis=0,**kwds):
  # x is an ndarray, with time along taxis
  x2 = np.moveaxis(x,taxis,0).reshape((len(t),-1))
  for i in range(x2.shape[1]):
    plt.plot(t,x2[:,i],**kwds)

def ribbon(t,x,taxis=0,interval=.9,alpha=.3,median=True,**kwds):
  # x is a list of ndarrays, with time along taxis
  x3 = np.stack([np.moveaxis(xi,taxis,0).reshape((len(t),-1)) for xi in x])
  qs = [interval_qs(i) for i in flatten(interval)]
  label = kwds.pop('label',None)
  for i in range(x3.shape[2]):
    for q in qs:
      plt.fill_between(t,
        np.nanquantile(x3[:,:,i],q[0],axis=0),
        np.nanquantile(x3[:,:,i],q[1],axis=0),
        **kwds,lw=.001,alpha=alpha)
    if median:
      plt.plot(t,np.nanquantile(x3[:,:,i],.5,axis=0),label=label,**kwds)

def boxplot(t,x,dt=5,tb=None,taxis=0,alpha=.3,width=.6,**kwds):
  # x is a list of ndarrays, with time along taxis
  if tb is None:
    tb = tdt(t,dt)
  else:
    dt = min(np.diff(dt))
  it = itslice(tb,t)
  g,Ng = kwds.pop('dodge',(0,1))
  tbd = np.array(tb) + dt * width * (g/Ng - .5 * (Ng-1)/Ng)
  x3 = np.stack([np.moveaxis(xi,taxis,0).reshape((len(t),-1)) for xi in x])
  label = kwds.pop('label',None)
  color = kwds.pop('color',(0,0,0))
  kwds.update(dict( # :\ really plt?
    medianprops  = dict(lw=.5,color=color),
    boxprops     = dict(lw=.5,ec=color,fc=(*color,alpha)),
    whiskerprops = dict(lw=.5,color=color),
    capprops     = dict(lw=.5,color=color),
    flierprops   = dict(lw=.5,mec=color,marker='+',mew=.5,ms=3),
  ))
  for i in range(x3.shape[2]):
    plt.plot(np.nan,np.nan,color=color,label=label) # for legend
    plt.boxplot(x3[:,it,i],positions=tbd,widths=width*.8*dt/Ng,patch_artist=True,**kwds)
  if dt < 10:
    plt.xticks(tdt(t,10),tdt(t,10))
  else:
    plt.xticks(tb,tb)

def ribbon_or_box(t,x,box=False,**kwds):
  if box:
    kwds.pop('ls',None)
    boxplot(t,x,dt=box,**kwds)
  else:
    kwds.pop('dodge',None)
    ribbon(t,x,**kwds)

def target(Ti,interval=.95,**kwds):
  label  = kwds.pop('label',None)
  ls     = kwds.pop('ls','-' if Ti.weight==1 else '--' if Ti.weight else ':')
  marker = kwds.pop('marker','D' if Ti.weight else 's')
  t = Ti.ind['t'] if Ti.vsop is None else Ti.ind1['t']
  Tm = Ti.mean()
  Terr = np.abs(np.reshape(Ti.ci(interval),(2,1)) - Tm)
  eh = plt.errorbar(t,Tm,yerr=Terr,marker=marker,ms=2,lw=1,mew=1,capsize=2,label=label,**kwds)
  eh[-1][0].set_linestyle(ls)

# ------------------------------------------------------------------------------
# plot function iterators

def targets_S(T,oname,skey,label=True,**kwds):
  if T is None: return
  S = strats[skey]
  for Ti in filter_targets(T,name=oname,ind=S.ind):
    target(Ti,color=S.color,label=(S.label if label else None),**kwds)
    label = False # only label once

def targets_vS(T,oname,skey1,skey2,vsop,label=True,**kwds):
  if T is None: return
  S1 = strats[skey1]
  S2 = strats[skey2]
  labelstr = out.vs_label(S1.label,S2.label,vsop)
  for Ti in filter_targets(T,name=oname,ind1=S1.ind,ind2=S2.ind,vsop=vsop):
    target(Ti,color=clr_interp(S1.color,S2.color),label=(labelstr if label else None),**kwds)
    label = False # only label once

def plot_S(fun,t,R,skey,box=False,**kwds):
  S = strats[skey]
  kwds.update(color=kwds.pop('color',S.color)) # TODO: make decorator?
  kwds.update(label=kwds.pop('label',S.label))
  fkwds = dict_split(kwds,['tvec','rate','t0']) # TODO: make decorator?
  if isinstance(fun,str):
    fun = out.by_name(fun)
  if isinstance(R,list):
    x = [fun(Ri,**S.ind,**fkwds) for Ri in R]
    ribbon_or_box(t,x,box=box,**kwds)
  else:
    x = fun(R,**S.ind,**fkwds)
    line(t,x,**kwds)

def plot_vS(fun,t,R,skey1,skey2,vsop,box=False,**kwds):
  S1 = strats[skey1]
  S2 = strats[skey2]
  kwds.update(color=kwds.pop('color',clr_interp(S1.color,S2.color)))
  kwds.update(label=kwds.pop('label',out.vs_label(S1.label,S2.label,vsop)))
  fkwds = dict_split(kwds,['tvec','rate','t0'])
  if isinstance(R,list):
    x = [out.vs_ind(fun,Ri,S1.ind,S2.ind,vsop,**fkwds) for Ri in R]
    ribbon_or_box(t,x,box=box,**kwds)
  else:
    x = out.vs_ind(fun,R,S1.ind,S2.ind,vsop,**fkwds)
    line(t,x,**kwds)

def plot_SvR(fun,t,R1,R2,skey,vsop,box=False,**kwds):
  S = strats[skey]
  kwds.update(color=kwds.pop('color',S.color))
  kwds.update(label=kwds.pop('label',S.label))
  fkwds = dict_split(kwds,['tvec','rate','t0'])
  if isinstance(R1,list):
    x = [out.vs_R(fun,R1i,R2i,vsop,**S.ind,**fkwds) for R1i,R2i in zip(R1,R2)]
    ribbon_or_box(t,x,box=box,**kwds)
  else:
    x = out.vs_R(fun,R1,R2,vsop,**S.ind,**fkwds)
    line(t,x,color=color,**kwds)
