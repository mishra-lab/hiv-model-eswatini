import numpy as np
from utils import fio,squarish,flatten
from model import system,target,out,plot,slicers

plotsize = 3 # inches
ttfname = fio.tmpfile('fit-{}.pdf')
specs = dict(
  NX            = dict(oname='NX',snames=['all']),
  Psi           = dict(oname='Psi',snames=['fsw.h','fsw.l','cli.h','cli.l']),
  Ph            = dict(oname='Ph',snames=['ahi','>500','<500','<350','<200'],ymax=1),
  prevalence    = dict(oname='prevalence',ymax=[.5,.5,.5,1]),
  prevalence1v2 = dict(oname='prevalence',snames=[('fsw.h','fsw.l'),('fsw','w'),('wh','wl'),('mh','ml')],vsop='1/2',ymax=5),
  incidence     = dict(oname='incidence',ymax=[.1,.1,.1,3]),
  incidence1v2  = dict(oname='incidence',snames=[('wh','wl'),('mh','ml')],vsop='1/2',ymax=100),
  diagnosed     = dict(oname='diagnosed',ymax=1),
  treated_c     = dict(oname='treated_c',ymax=1),
  treated_u     = dict(oname='treated_u',ymax=1),
  vls_c         = dict(oname='vls_c',ymax=1),
  vls_u         = dict(oname='vls_u',ymax=1),
  condom        = dict(oname='condom',snames=['msp','cas','swo','swr'],ymax=1),
  circum        = dict(oname='circum',snames=['*']),
  dx_rate       = dict(oname='dx_rate'),
  tx_rate       = dict(oname='tx_rate'),
)
specsets = dict(
  hiv     = ['prevalence','prevalence1v2','incidence','incidence1v2'],
  cascade = ['diagnosed','treated_c','treated_u','vls_c','vls_u'],
  extra   = ['NX','Psi','Ph','condom','circum'],
  rates   = ['dx_rate','tx_rate'],
)

def plot_sets(t,Rs,T=None,fname='pyplots.pdf',debug=False,sets=None,snames=None):
  if sets is None: sets = specsets.keys()
  if snames is None: snames = ['all','w','m','fsw']
  kwds = dict(T=T,snames=snames)
  if debug: # best 100%, 10%, 1% fits as ribbons; then merge 
    Rs = system.drop_fails(Rs)[0]
    Rss = [target.top_ll(Rs,top) for top in (1.,.1,.01)]
  else: # 100% fits as ribbon + median; no merge
    Rss = [Rs]
    kwds.update(tfname=fname)
  tfnames = [plot_output(t,Rss,**dict(kwds,**specs[name])) \
    for set in flatten(sets) for name in specsets[set]]
  if debug: fio.pdfmerge(fname,tfnames)

def plot_output(t,Rss,oname,snames,T=None,tfname=None,ylab=None,ymax=None,**kwds):
  if tfname is None: tfname = ttfname
  if ylab is None: ylab = out.labels.get(oname,oname)
  if np.size(ymax) == 1: ymax = len(snames) * flatten(ymax)
  fh,ah = plot.subplots(1,len(snames))
  kwds.update(interval=1,median=(len(Rss)==1))
  for s,sname in enumerate(snames):
    plot.plt.sca(ah[0,s])
    if isinstance(sname,tuple): # 2-group (vs) type
      title = out.vs_label(slicers[sname[0]].label,slicers[sname[1]].label,'\n')
      plot.labels(title=title,x=None,y=ylab+' Ratio' if s==0 else None)
      for Rs in Rss: # ribbons
        plot.plot_vS(oname,t,Rs,sname[0],sname[1],**kwds)
      plot.targets_vS(T,oname,sname[0],sname[1],kwds['vsop'])
    else: # 1-group type
      plot.labels(title=slicers[sname].label,x=None,y=ylab if s==0 else None)
      for Rs in Rss: # ribbons
        plot.plot_S(oname,t,Rs,sname,**kwds)
      plot.targets_S(T,oname,sname)
    plot.plt.ylim((0,ymax[s]))
  fh.set_size_inches((.5+plotsize*len(snames),plotsize))
  fh.tight_layout()
  return plot.save(tfname.format(oname+(kwds.get('vsop','').replace('/','v'))))
