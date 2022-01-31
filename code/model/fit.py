import numpy as np
from utils import flatten,squarish,fio
from model import slicers,system,out,target,plot

plotsize = 3 # inches
tfname = '.fit-{}.pdf'

def plot_all(t,Rs,T,fname='fit.pdf',tops=(.5,.10,.02),drop=True):
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
    plot_output(t,Rss,'NX',        ['ALL'],T=T),
    plot_output(t,Rss,'Ph',        ['AHI','>500','<500','<350','<200'],T=T),
    plot_output(t,Rss,'X_rate',    ['ALL'],rate='death_hc',ylab='HIV Mortality'),
    plot_output(t,Rss,'prevalence',['ALL','W','M','FSW'],T=T,ylim=(0,1)),
    plot_output(t,Rss,'prevalence',
      [('FSW.H','FSW.L'),('FSW','W'),('WH','WL'),('Cli.H','Cli.L'),('Cli','M'),('MH','ML')],
      vsop='1/2',T=T,ylim=(1,5)),
    plot_output(t,Rss,'incidence', ['ALL','W','M','FSW'],T=T),
    plot_output(t,Rss,'diagnosed', ['ALL','W','M','FSW'],T=T),
    plot_output(t,Rss,'treated_c', ['ALL','W','M','FSW'],T=T),
    plot_output(t,Rss,'vls_c',     ['ALL','W','M','FSW'],T=T),
    plot_output(t,Rss,'treated_u', ['ALL','W','M','FSW'],T=T),
    plot_output(t,Rss,'vls_u',     ['ALL','W','M','FSW'],T=T),
    plot_output(t,Rss,'condom',    ['LT','ST','SWR','SWO']),
    plot_output(t,Rss,'circum',    ['*']),
    plot_output(t,Rss,'gud',       ['*']),
    # plot_output(t,Rss,'dx_rate',   ['W','M','FSW']),
    # plot_output(t,Rss,'tx_rate',   ['W','M','FSW']),
    # plot_output(t,Rss,'tx_rate',   ['W','M','FSW']),
    # plot_output(t,Rss,'cuminfect', ['W','M','FSW'],tvec=t,T=T),
  ]
  fio.pdfmerge(fname,tfnames,rm=True)

def plot_param(Rss,pname,t=None,**kwds):
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

def plot_output(t,Rss,oname,snames,T=None,ylab=None,ylim=None,**kwds):
  if ylab is None: ylab = out.labels.get(oname,oname)
  fh,ah = plot.subplots(1,len(snames))
  kwds.update(median=False,interval=1)
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
