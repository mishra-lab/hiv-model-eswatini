# functions for plotting model outputs vs calibration targets

import numpy as np
from utils import fio,flatten
from model import system,target,out,plot,strats

plotsize = 3 # inches
specs = dict(
  NX      = dict(oname='NX',skeys=['all']),
  Psi     = dict(oname='Psi',skeys=['fsw.h','fsw.l','cli.h','cli.l']),
  Ph      = dict(oname='Ph',skeys=['ahi','>500','<500','<350','<200'],ymax=1),
  prev    = dict(oname='prevalence',ymax=[.5,.5,.5,1]),
  prev1v2 = dict(oname='prevalence',skeys=[('fsw.h','fsw.l'),('fsw','w'),('wh','wl'),('mh','ml')],vsop='1/2',ymax=5),
  prevanc = dict(oname='prevalence',skeys=['w'],ymax=[.5],T=target.get_prevalence_esw_anc()),
  inc     = dict(oname='incidence',ymax=[.1,.1,.1,2]),
  inc1v2  = dict(oname='incidence',skeys=[('wh','wl'),('mh','ml')],vsop='1/2',ymax=100),
  cuminf  = dict(oname='cuminfect'),
  tdsc    = dict(oname='tdsc',skeys=['msp','cas','swo','swr'],ymax=1),
  diag    = dict(oname='diagnosed',ymax=1),
  treat_c = dict(oname='treated_c',ymax=1),
  treat_u = dict(oname='treated_u',ymax=1),
  vls_c   = dict(oname='vls_c',ymax=1),
  vls_u   = dict(oname='vls_u',ymax=1),
  condom  = dict(oname='condom',skeys=['msp','cas','swo','swr'],ymax=1),
  circum  = dict(oname='circum',skeys=['*']),
  dx_rate = dict(oname='dx_rate'),
  tx_rate = dict(oname='tx_rate'),
)
specsets = dict(
  hiv     = ['prev','prev1v2','inc','inc1v2'],
  cascade = ['diag','treat_c','treat_u','vls_c','vls_u'],
  extra   = ['NX','Psi','Ph','tdsc','condom','circum','prevanc'],
  rates   = ['dx_rate','tx_rate'],
)

def plot_sets(t,Rs,T=None,tfname=None,debug=False,sets=None,skeys=None):
  # helper for plotting groups of figures
  if tfname is None: tfname = fio.tmpfile('fit-{}.pdf')
  if sets is None: sets = specsets.keys()
  if skeys is None: skeys = ['all','w','m','fsw']
  kwds = dict(T=T,skeys=skeys)
  if debug: # best 100%, 10%, 1% fits as ribbons
    Rss = [target.top_ll(Rs,top) for top in (1.,.1,.01)]
  else: # 100% fits as ribbon + median
    Rss = [Rs]
  tfnames = [plot_output(t,Rss,**dict(kwds,fname=tfname.format(name),**specs[name]))
    for set in flatten(sets) for name in specsets[set]
    if name != 'tdsc' or out.can_tdsc(Rs[0])] # can't plot tdsc if didn't save Xk
  if debug: fio.pdfmerge('pyplots.pdf',tfnames) # merge figs into 1 pdf

def plot_output(t,Rss,oname,skeys,fname,T=None,ylab=None,ymax=None,**kwds):
  # plot output (ribbons) for multiple strata, plus any matching targets
  if ylab is None: ylab = out.labels.get(oname,oname)
  if np.size(ymax) == 1: ymax = len(skeys) * flatten(ymax)
  fh,ah = plot.subplots(1,len(skeys)) # subplot per stratum
  kwds.update(interval=1,median=(len(Rss)==1)) # HACK: plot median if not debug
  if oname == 'cuminfect': kwds.update(tvec=t) # need original t vec for cuminf
  for s,skey in enumerate(skeys): # strata
    plot.plt.sca(ah[0,s])
    if isinstance(skey,tuple): # 2-group (vs) type
      title = out.vs_label(strats[skey[0]].label,strats[skey[1]].label,'\n')
      plot.labels(title=title,x=None,y=ylab+' Ratio' if s==0 else None)
      for Rs in Rss: # ribbons
        plot.plot_vS(oname,t,Rs,skey[0],skey[1],**kwds)
      plot.targets_vS(T,oname,skey[0],skey[1],kwds['vsop'])
    else: # 1-group type
      plot.labels(title=strats[skey].label,x=None,y=ylab if s==0 else None)
      for Rs in Rss: # ribbons
        plot.plot_S(oname,t,Rs,skey,**kwds)
      plot.targets_S(T,oname,skey)
    plot.plt.ylim((0,ymax[s]))
  fh.set_size_inches((.5+plotsize*len(skeys),plotsize))
  fh.tight_layout()
  return plot.save(fname)
