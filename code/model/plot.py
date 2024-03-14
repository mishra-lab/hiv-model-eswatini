# functions for plotting model outputs & calibration targets

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from model.target import filter_targets
from utils import _,log,genpath,flatten,dict_split,interval_qs,clr_interp,itslice,tdt
from model import out,strats

# ------------------------------------------------------------------------------
# helper functions
# TODO: some of these should go in utils (?)

def subplots(row,col):
  # initialize a grid of subplots & return handles
  fh,ah = plt.subplots(row,col) # figure & axes handles
  ah = np.reshape(ah,(row,col))
  return fh,ah

def save(fname,*args,close=True,**kwds):
  # save a plot with some helpful default actions
  log(2,'plot.save: '+fname)
  plt.savefig(genpath(fname),*args,**kwds)
  if close: plt.close()
  return fname

def labels(x='Year',y=None,title=None,fs=11):
  # conveniently label x,y,title
  if x     is not None: plt.xlabel(x,fontsize=fs)
  if y     is not None: plt.ylabel(y,fontsize=fs)
  if title is not None: plt.title(title,fontsize=fs)

# ------------------------------------------------------------------------------
# main plot functions

def line(t,x,**kwds):
  # simple line plot: x vs t (vectors)
  plt.plot(t,x,**kwds)

def ribbon(t,xs,interval=.9,alpha=.3,median=True,**kwds):
  # plot the range of xs vs t, where xs is a list of vectors
  x2 = np.stack(xs)
  qs = [interval_qs(i) for i in flatten(interval)]
  label = kwds.pop('label',None)
  for q in qs:
    plt.fill_between(t,
      np.nanquantile(x2,q[0],axis=0),
      np.nanquantile(x2,q[1],axis=0),
      **kwds,lw=.001,alpha=alpha)
  if median:
    line(t,np.nanquantile(x2,.5,axis=0),label=label,**kwds)

def boxplot(t,xs,dt=5,tb=None,alpha=.3,width=.6,**kwds):
  # boxplots of xs vs t, at regular intervals (dt) or at specific times (tb)
  # xs is a list of vectors as above
  if tb is None: tb = tdt(t,dt)  # define tb from dt
  else: dt = min(np.diff(tb))    # define dt from tb (for spacing only)
  it = itslice(tb,t)             # get indices of tb in t
  g,Ng = kwds.pop('dodge',(0,1)) # TODO
  tbp = np.array(tb) + dt * width * (g/Ng - .5 * (Ng-1)/Ng) # x positions
  x2 = np.stack(xs)
  label = kwds.pop('label',None)
  color = kwds.pop('color',(0,0,0))
  kwds.update(dict( # define style, no shortcut apparently
    medianprops  = dict(lw=.5,color=color),
    boxprops     = dict(lw=.5,ec=color,fc=(*color,alpha)),
    whiskerprops = dict(lw=.5,color=color),
    capprops     = dict(lw=.5,color=color),
    flierprops   = dict(lw=.5,mec=color,marker='+',mew=.5,ms=3),
  ))
  plt.plot(np.nan,np.nan,color=color,label=label) # HACK for legend
  plt.boxplot(x2[:,it],positions=tbp,widths=width*.8*dt/Ng,patch_artist=True,**kwds)
  if dt < 10: plt.xticks(tdt(t,10),tdt(t,10))
  else:       plt.xticks(tb,tb)

def ribbon_or_box(t,xs,box=False,**kwds):
  # helper to plot xs vs t as ribbon (default) or boxplot
  if box:
    kwds.pop('ls',None) # remove unused 'ls' kwds
    boxplot(t,xs,dt=box,**kwds)
  else:
    kwds.pop('dodge',None) # remove unused 'dodge' kwds
    ribbon(t,xs,**kwds)

def target(Ti,interval=.95,**kwds):
  # plot Target Ti mean & 95% CI as point & errorbar
  label  = kwds.pop('label',None)
  ls     = kwds.pop('ls','-' if Ti.weight==1 else '--' if Ti.weight else ':')
  marker = kwds.pop('marker','D' if Ti.weight else 's') # diamond if weight > 0 else square
  t = Ti.ind['t'] if Ti.vsop is None else Ti.ind1['t']  # target time
  Tm = Ti.mean()
  Terr = np.abs(np.reshape(Ti.ci(interval),(2,1)) - Tm) # 95% CI relative to Tm for below
  eh = plt.errorbar(t,Tm,yerr=Terr,marker=marker,ms=2,lw=1,mew=1,capsize=2,label=label,**kwds)
  eh[-1][0].set_linestyle(ls) # HACK to specify linestyle (ls)

# ------------------------------------------------------------------------------
# plot function iterators

def targets_S(T,oname,skey,label=True,**kwds):
  # plot corresponding targets for this oname & skey
  if T is None: return
  S = strats[skey]
  for Ti in filter_targets(T,name=oname,ind=S.ind):
    target(Ti,color=S.color,label=(S.label if label else None),**kwds)
    label = False # only label once

def targets_vS(T,oname,skey1,skey2,vsop,label=True,**kwds):
  # plot corresponding vs-type targets for this oname & skeys
  if T is None: return
  S1 = strats[skey1]
  S2 = strats[skey2]
  labelstr = out.vs_label(S1.label,S2.label,vsop)
  for Ti in filter_targets(T,name=oname,ind1=S1.ind,ind2=S2.ind,vsop=vsop):
    target(Ti,color=clr_interp(S1.color,S2.color),label=(labelstr if label else None),**kwds)
    label = False # only label once

def plot_S(fun,t,R,skey,box=False,**kwds):
  # plot an output for this skey and 1+ model results (R)
  # fun can be a string like oname or the output-computing function
  S = strats[skey]
  kwds.update(color=kwds.pop('color',S.color))
  kwds.update(label=kwds.pop('label',S.label))
  fkwds = dict_split(kwds,['tvec','rate','t0']) # pop non-plotting kwds
  if isinstance(fun,str): fun = out.by_name(fun)
  if isinstance(R,list):
    xs = [fun(Ri,**S.ind,**fkwds) for Ri in R]
    ribbon_or_box(t,xs,box=box,**kwds)
  else:
    x = fun(R,**S.ind,**fkwds)
    line(t,x,**kwds)

def plot_vS(fun,t,R,skey1,skey2,vsop,box=False,**kwds):
  # plot a vs-type output for these skeys and 1+ model results (R)
  # fun can be a string like oname or the output-computing function
  S1 = strats[skey1]
  S2 = strats[skey2]
  kwds.update(color=kwds.pop('color',clr_interp(S1.color,S2.color)))
  kwds.update(label=kwds.pop('label',out.vs_label(S1.label,S2.label,vsop)))
  fkwds = dict_split(kwds,['tvec','rate','t0'])
  if isinstance(R,list):
    xs = [out.vs_ind(fun,Ri,S1.ind,S2.ind,vsop,**fkwds) for Ri in R]
    ribbon_or_box(t,xs,box=box,**kwds)
  else:
    x = out.vs_ind(fun,R,S1.ind,S2.ind,vsop,**fkwds)
    line(t,x,**kwds)

def plot_SvR(fun,t,R1,R2,skey,vsop,box=False,**kwds):
  # vs-plot an output for this skey and pair(s) of 1+ model results R1 vs R2
  # fun can be a string like oname or the output-computing function
  S = strats[skey]
  kwds.update(color=kwds.pop('color',S.color))
  kwds.update(label=kwds.pop('label',S.label))
  fkwds = dict_split(kwds,['tvec','rate','t0'])
  if isinstance(R1,list):
    xs = [out.vs_R(fun,R1i,R2i,vsop,**S.ind,**fkwds) for R1i,R2i in zip(R1,R2)]
    ribbon_or_box(t,xs,box=box,**kwds)
  else:
    x = out.vs_R(fun,R1,R2,vsop,**S.ind,**fkwds)
    line(t,x,color=color,**kwds)
