import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from model.target import filter_targets
from utils import genpath,flatten,clr_interp,interval_qs,squarish,dict_split,itslice,tdt,log,globs
from model import _,dimkeys,dimensions,slicers,out

# TODO: labeling & colors for vsop

# HELPERS ----------------------------------------------------------------------
# TODO: some of these should go in utils

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

def sliceiter(shape,dstr,ah=None,join='\n'):
  # iterate over all stratifications for a shape;
  # and given dimension "letters" in dstr, build up the associated names (from model/__init__.py)
  # e.g. sliceiter([2,4],'si') -> 3rd yield = (2,(0,2),'Women\nLower Risk Sex Work')
  # also: initialize (or reuse) an array of subplots, and plt.sca() before each yield;
  # and if we have too many subplots in the grid, turn unused axes off
  size   = int(np.prod(shape))
  grid   = np.moveaxis(np.indices(shape),0,-1)
  slices = [tuple(i) for i in grid.reshape(size,-1)]
  if dstr is None:
    dlabfun = lambda i: str(slices[i])
  else:
    dkeys = list(dimensions.keys())
    dlabs = [dimensions[dkeys[dimkeys.index(d)]] for d in dstr]
    dlabfun = lambda i: join.join([dlabs[d][slices[i][d]] for d in range(len(dstr))])
  if ah is None:
    ah = subplots(*squarish(size))[1]
  for i,ahi in enumerate(flatten(ah)):
    if i < size:
      plt.sca(ahi)
      yield i,slices[i],dlabfun(i)
    else:
      ahi.set_axis_off()

# MAIN PLOT FUNCTIONS ----------------------------------------------------------

def line(t,x,taxis=0,**kwds):
  # x is an ndarray, with time along taxis
  x2 = np.moveaxis(x,taxis,0).reshape((len(t),-1))
  for i in range(x2.shape[1]):
    plt.plot(t,x2[:,i],**kwds)

def ribbon(t,x,taxis=0,interval=.9,alpha=.2,median=True,**kwds):
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

def boxplot(t,x,dt=5,tb=None,taxis=0,alpha=.2,width=.6,**kwds):
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
  t = Ti.pop['t'] if Ti.vsop is None else Ti.pop1['t']
  Tm = Ti.mean()
  Terr = np.abs(np.reshape(Ti.ci(interval),(2,1)) - Tm)
  eh = plt.errorbar(t,Tm,yerr=Terr,marker=marker,ms=2,lw=1,mew=1,capsize=2,label=label,**kwds)
  eh[-1][0].set_linestyle(ls)

def hist_p(Ps,name,dstr=None,ah=None,alpha=.3,xlim=None,gbins=False,**kwds):
  def ihist(x,title='',**kwds):
    bins = plt.hist(x,alpha=alpha,**kwds)[1]
    plt.title(title,loc='left',y=1,pad=-12*(1+title.count('\n')))
    plt.yticks([])
    plt.xlim(xlim)
    return bins
  gbkey = 'hist_p_'+name+'_{}'
  Xs = [np.squeeze(P[name]) for P in Ps]
  shape = np.shape(Xs[0])
  if len(shape):
    for i,slicer,ilabel in sliceiter(shape,dstr,ah):
      if gbins: kwds['bins'] = globs(get=gbkey.format(i),default=32)
      bins = ihist([X[slicer] for X in Xs],ilabel,**kwds)
      if gbins: globs(**{gbkey.format(i):bins})
  else:
    if ah is None: ah = subplots(1,1)[1]
    if gbins: kwds['bins'] = globs(get=gbkey.format(0),default=32)
    bins = ihist(Xs,**kwds)
    if gbins: globs(**{gbkey.format(0):bins})
  plt.suptitle(name)
  return ah

# PLOT FUNCTION ITERATORS ------------------------------------------------------

def targets_S(T,oname,sname,label=True,**kwds):
  if T is None: return
  S = slicers[sname]
  for Ti in filter_targets(T,name=oname,pop=S.pop):
    target(Ti,color=S.color,label=(S.label if label else None),**kwds)
    label = False # only label once

def targets_vS(T,oname,sname1,sname2,vsop,label=True,**kwds):
  if T is None: return
  S1 = slicers[sname1]
  S2 = slicers[sname2]
  labelstr = out.vs_label(S1.label,S2.label,vsop)
  for Ti in filter_targets(T,name=oname,pop1=S1.pop,pop2=S2.pop,vsop=vsop):
    target(Ti,color=clr_interp(S1.color,S2.color),label=(labelstr if label else None),**kwds)
    label = False # only label once

def plot_S(fun,t,R,sname,box=False,**kwds):
  S = slicers[sname]
  kwds.update(color=kwds.pop('color',S.color)) # TODO: make decorator?
  kwds.update(label=kwds.pop('label',S.label))
  fkwds = dict_split(kwds,['tvec','rate','t0']) # TODO: make decorator?
  if isinstance(fun,str):
    fun = out.by_name(fun)
  if isinstance(R,list):
    x = [fun(Ri,**S.pop) for Ri in R]
    ribbon_or_box(t,x,box=box,**kwds)
  else:
    x = fun(R,**S.pop,**fkwds)
    line(t,x,**kwds)

def plot_vS(fun,t,R,sname1,sname2,vsop,box=False,**kwds):
  S1 = slicers[sname1]
  S2 = slicers[sname2]
  kwds.update(color=kwds.pop('color',clr_interp(S1.color,S2.color)))
  kwds.update(label=kwds.pop('label',out.vs_label(S1.label,S2.label,vsop)))
  fkwds = dict_split(kwds,['tvec','rate','t0'])
  if isinstance(R,list):
    x = [out.vs_pop(fun,Ri,S1.pop,S2.pop,vsop,**fkwds) for Ri in R]
    ribbon_or_box(t,x,box=box,**kwds)
  else:
    x = out.vs_pop(fun,R,S1.pop,S2.pop,vsop,**fkwds)
    line(t,x,**kwds)

def plot_SvR(fun,t,R1,R2,sname,vsop,box=False,**kwds):
  S = slicers[sname]
  kwds.update(color=kwds.pop('color',S.color))
  kwds.update(label=kwds.pop('label',S.label))
  fkwds = dict_split(kwds,['tvec','rate','t0'])
  if isinstance(R1,list):
    x = [out.vs_R(fun,R1i,R2i,vsop,**S.pop,**fkwds) for R1i,R2i in zip(R1,R2)]
    ribbon_or_box(t,x,box=box,**kwds)
  else:
    x = out.vs_R(fun,R1,R2,vsop,**S.pop,**fkwds)
    line(t,x,color=color,**kwds)

# TODO: labels for vsop
