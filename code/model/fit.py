import numpy as np
from utils import flatten,squarish,fio
from model import slicers,system,out,target,plot

plotsize = 3 # inches
ttfname = '.fit-{}.pdf'

def plot_all(t,Rs,T,fname='fit.pdf',tops=(1.,.2,.04),drop=True):
  if drop: Rs = system.drop_fails(Rs)[0]
  Rss = [target.top_q_ll(Rs,top) for top in tops] if tops else [Rs]
  tfnames = [
    # param histograms
    plot_param(Rss,'PX_si',    dstr='si'),
    plot_param(Rss,'A_ap',     dstr='ap'),
    plot_param(Rss,'C_psi',    dstr='psi'),
    plot_param(Rss,'Rbeta_as', dstr='as'),
    plot_param(Rss,'Rbeta_h',  dstr='h'),
    plot_param(Rss,'P_gud',    dstr='si'),
    plot_param(Rss,'EHY_acute',dstr=''),
    plot_param(Rss,'t0_hiv',   dstr=''),
    # output projections
    plot_output(t,Rss,'NX',        ['all'],T=T),
    plot_output(t,Rss,'Ph',        ['ahi','>500','<500','<350','<200'],T=T),
    plot_output(t,Rss,'X_rate',    ['all'],rate='death_hc',ylab='HIV Mortality',T=T),
    plot_output(t,Rss,'prevalence',['all','w','m','fsw'],T=T,ylim=(0,1)),
    plot_output(t,Rss,'prevalence',
      [('fsw.h','fsw.l'),('fsw','w'),('wh','wl'),('cli.h','cli.l'),('cli','m'),('mh','ml')],
      vsop='1/2',T=T,ylim=(1,5)),
    plot_output(t,Rss,'incidence', ['all','w','m','fsw'],T=T),
    plot_output(t,Rss,'diagnosed', ['all','w','m','fsw'],T=T),
    plot_output(t,Rss,'treated_c', ['all','w','m','fsw'],T=T),
    plot_output(t,Rss,'vls_c',     ['all','w','m','fsw'],T=T),
    plot_output(t,Rss,'treated_u', ['all','w','m','fsw'],T=T),
    plot_output(t,Rss,'vls_u',     ['all','w','m','fsw'],T=T),
    plot_output(t,Rss,'condom',    ['msp','cas','swo','swr'],T=T),
    plot_output(t,Rss,'circum',    ['*'],T=T),
    plot_output(t,Rss,'gud',       ['*'],T=T),
    # plot_output(t,Rss,'dx_rate',   ['W','M','FSW']),
    # plot_output(t,Rss,'tx_rate',   ['W','M','FSW']),
    # plot_output(t,Rss,'tx_rate',   ['W','M','FSW']),
    # plot_output(t,Rss,'cuminfect', ['W','M','FSW'],tvec=t,T=T),
  ]
  fio.pdfmerge(fname,tfnames,rm=True)

def plot_fit(t,Rs,T,fname,tops=(1.,.2,.04),drop=True,merge=True):
  if drop: Rs = system.drop_fails(Rs)[0]
  Rss = [target.top_q_ll(Rs,top) for top in tops] if (tops and T) else [Rs]
  kwds = dict(T=T,tfname=(None if merge else fname))
  tfnames = [
    plot_output(t,Rss,'prevalence',['all','w','m','fsw'],**kwds,ylim=(0,1)),
    plot_output(t,Rss,'incidence', ['all','w','m','fsw'],**kwds),
    plot_output(t,Rss,'diagnosed', ['all','w','m','fsw'],**kwds,ylim=(0,1)),
    plot_output(t,Rss,'treated_c', ['all','w','m','fsw'],**kwds,ylim=(0,1)),
    plot_output(t,Rss,'vls_c',     ['all','w','m','fsw'],**kwds,ylim=(0,1)),
    plot_output(t,Rss,'condom',    ['msp','cas','swo','swr'],**kwds,ylim=(0,1)),
    plot_output(t,Rss,'circum',    ['*'],**kwds,ylim=(0,1)),
    plot_output(t,Rss,'NX',        ['all'],**kwds),
    plot_output(t,Rss,'Ph',        ['ahi','>500','<500','<350','<200'],**kwds),
    plot_output(t,Rss,'prevalence',
      [('fsw.h','fsw.l'),('fsw','w'),('cli.h','cli.l'),('cli','m')],
      vsop='1/2',**kwds,ylim=(1,5)),
  ]
  if merge:
    fio.pdfmerge(fname,tfnames,rm=True)

def plot_refit(t,Rs,T,fname,tops=(1.,.2,.04),drop=True,merge=True):
  if drop: Rs = system.drop_fails(Rs)[0]
  Rss = [target.top_q_ll(Rs,top) for top in tops] if (tops and T) else [Rs]
  kwds = dict(T=T,tfname=(None if merge else fname))
  groups = ['all','aq','fsw','cli']
  tfnames = [
    # output projections
    plot_output(t,Rss,'prevalence',groups,**kwds,ylim=(0,1)),
    plot_output(t,Rss,'incidence', groups,**kwds),
    plot_output(t,Rss,'diagnosed', groups,**kwds,ylim=(0,1)),
    plot_output(t,Rss,'treated_c', groups,**kwds,ylim=(0,1)),
    plot_output(t,Rss,'vls_c',     groups,**kwds,ylim=(0,1)),
    plot_output(t,Rss,'treated_u', groups,**kwds,ylim=(0,1)),
    plot_output(t,Rss,'vls_u',     groups,**kwds,ylim=(0,1)),
  ]
  if merge:
    fio.pdfmerge(fname,tfnames,rm=True)

def plot_param(Rss,pname,tfname=None,**kwds):
  if tfname is None: tfname = ttfname
  p0 = Rss[0][0]['P'][pname]
  row,col = squarish(np.size(p0))
  fh,ah = plot.subplots(row,col)
  kwds.update(ah=ah,gbins=True,density=True)
  cmap = plot.cmap(len(Rss),trim=True)
  for r,Rs in enumerate(Rss): # subsets of model fits
    plot.hist_p([R['P'] for R in Rs],pname,color=cmap[r],**kwds)
  fh.set_size_inches((plotsize*col,plotsize*row))
  fh.tight_layout()
  return plot.save(tfname.format(pname))

def plot_output(t,Rss,oname,snames,T=None,tfname=None,ylab=None,ylim=None,**kwds):
  if tfname is None: tfname = ttfname
  if ylab is None: ylab = out.labels.get(oname,oname)
  fh,ah = plot.subplots(1,len(snames))
  kwds.update(median=False,interval=1 if T else (1,.5,.1))
  for s,sname in enumerate(snames): # subplots
    plot.plt.sca(ah[0,s])
    if isinstance(sname,tuple):
      plot.labels(title=out.vs_label(slicers[sname[0]].label,slicers[sname[1]].label,kwds['vsop']),
        x=None,y=ylab if s==0 else None)
      for Rs in Rss:
        plot.plot_vS(oname,t,Rs,sname[0],sname[1],**kwds)
      plot.targets_vS(T,oname,sname[0],sname[1],kwds['vsop'])
    else:
      plot.labels(title=slicers[sname].label,x=None,y=ylab if s==0 else None)
      for Rs in Rss:
        plot.plot_S(oname,t,Rs,sname,**kwds)
      plot.targets_S(T,oname,sname)
    plot.plt.ylim(ylim)
  fh.set_size_inches((plotsize*len(snames),plotsize))
  fh.tight_layout()
  return plot.save(tfname.format(oname+(kwds.get('vsop','').replace('/','v'))))
