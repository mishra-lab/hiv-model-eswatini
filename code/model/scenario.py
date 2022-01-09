import numpy as np
from utils import tarray as ta
from utils import os,stats,rootpath
from utils.ops import flatten,ppool,minimize,npsave,npload,filehash
from model import slicers,params,system,target,out,plot,handfit
import ipdb

plotsize = 3 # inches
hash = filehash('params.py','system.py','target.py',root=rootpath('code','model'))

def Rqx_by_group(Rqx=None,**kwds):
  # create/update a relative-rate using pop=value kwds, e.g. FSW=1.5
  if Rqx is None: Rqx = np.ones((2,4,1))
  for sname,sRqx in kwds.items():
    for sk in flatten(slicers[sname].pop.get('s',(0,1))):
      for ik in flatten(slicers[sname].pop.get('i',(0,1,2,3))):
        Rqx[sk,ik,0] = sRqx
  return Rqx

def fit_cascade(P,PD,T,t,P0=None,ftol=.01):
  # fit the model to "T", with base parameters P, and re-sampled values in PD
  # PD has keys like "rate:pop"; we re-sample based on percentiles (pfit)
  if P0 is None: P0 = [.1 for i in PD]
  def jfun(pfits): # optimization function
    # map pfits (quantiles) -> values from PD & update P with the new values
    for pfit,key in zip(pfits,PD.keys()):
      rate,pop = key.split(':')
      P.update({rate:Rqx_by_group(Rqx=P.get(rate),**{pop:PD[key].ppf(pfit)})})
    # run the model
    R = system.run(P,t,T,RPts=[])
    # print((R['ll'],pfits),flush=True) # DEBUG
    return -R['ll']
  # run the optimization:
  # SLSQP works best; all bounds (0,1) due to percentile trick; ftol controls precision
  M = minimize(jfun,P0,method='SLSQP',bounds=len(PD)*((0,1),),options=dict(ftol=ftol))
  if M.success:
    return P # TODO: is this robust updated P?
  return False

def fit_cascade_n(Ps,PD,T,t,P0=None,**kwds):
  fun = lambda P: fit_cascade(P,PD=PD,T=T,t=t,P0=P0,**kwds)
  return ppool(min(len(Ps),7)).map(fun,Ps)

def plot_diff(t,R1s,R2s,oname,snames,intervals=.95,**kwds):
  fh,ah = plot.subplots(1,len(snames))
  intervals = flatten(intervals)
  for i,sname in enumerate(snames): # subplots
    plot.plt.sca(ah[0,i])
    plot.labels(title=slicers[sname].label,
      x='Time (years)',y=out.labels.get(oname,oname) if i==0 else None)
    for interval in intervals:
      plot.plot_SvR(oname,t,R1s,R2s,sname,vsop='1-2/1',tvec=t,interval=interval)
  fh.set_size_inches((plotsize*len(snames),plotsize))
  fh.tight_layout()
  return fh,ah

def fname(key,N,N0=0,ext=''):
  return '{}_N={}-{}{}'.format(key,N0,N0+N-1,ext)

def run(N,N0=0,sample=True,refit=True,top=.10):
  t  = system.f_t(t1=2025)
  t2 = system.f_t(t1=2050)
  T1 = target.get_all_esw()
  if sample:
    # run model with random samples
    P1s_sam = params.get_n_all(N,seeds=range(N0,N0+N))
    R1s_sam = system.run_n(P1s_sam,t,T1)
    R1s_sam = system.drop_fails(R1s_sam)[0]
    handfit.plot_all(t,R1s_sam,T1,fname=fname('handfit_sam',N,N0,'.pdf'))
    # save top % samples
    R1s = target.top_q_ll(R1s_sam,top=top)
    P1s = [R1['P'] for R1 in R1s]
    npsave(os.path.join(hash,fname('P1s',N,N0)),P1s)
  else:
    P1s = npload(os.path.join(hash,fname('P1s',N,N0)))
  R1s = system.run_n(P1s,t2,T1) # re-run for t2
  R1s = system.drop_fails(R1s)[0] # TODO: better solution?
  P1s = [R1['P'] for R1 in R1s]   # TODO: better solution?
  handfit.plot_all(t2,R1s,T1,fname=fname('handfit_base',N,N0,'.pdf'))
  # re-fitting definitions
  PDD = {'0-1':stats.uniform(0,1),'1-1':stats.uniform(1-1e-9,1+1e-9),'1-50':stats.uniform(1,50)}
  for case in ['LoLo','LoHi','HiLo']:
    if case=='LoLo':
      T2 = target.get_cascade_2020('ssa-lowest')+target.get_cascade_2020('ssa-lowest-fsw')
      PD = {'Rdx:ALL':PDD['0-1'],'Rtx:ALL':PDD['0-1'],'Rux:ALL':PDD['1-50'],
            'Rdx:FSW':PDD['0-1'],'Rtx:FSW':PDD['0-1'] }
    if case=='LoHi':
      T2 = target.get_cascade_2020('ssa-lowest')
      PD = {'Rdx:ALL':PDD['0-1'],'Rtx:ALL':PDD['0-1'],'Rux:ALL':PDD['1-50'],
            'Rdx:FSW':PDD['1-1'],'Rtx:FSW':PDD['1-1'] }
    if case=='HiLo':
      T2 = target.get_cascade_2020('ssa-lowest-fsw')
      PD = {'Rdx:ALL':PDD['1-1'],'Rtx:ALL':PDD['1-1'],'Rux:ALL':PDD['1-1'],
            'Rdx:FSW':PDD['0-1'],'Rtx:FSW':PDD['0-1'],'Rux:FSW':PDD['1-50'] }
    if refit:
      # resample relative-dx,tx,ux to fit T2
      P2s = fit_cascade_n(P1s,PD,T2,t,ftol=.01)
      npsave(os.path.join(hash,fname('P2s_'+case,N,N0)),P2s)
    else:
      P2s = npload(os.path.join(hash,fname('P2s_'+case,N,N0)))
    R2s = system.run_n(P2s,t2,T2) # re-run for t2
    handfit.plot_all(t2,R2s,T2,fname=fname('handfit_'+case,N,N0,'.pdf'))
    # plot differences
    fh,ah = plot_diff(t2,R2s,R1s,'cuminfect',['ALL','FSW','Cli'],intervals=[.5,.95])
    for ahi in ah.flatten(): ahi.set_ylim((0,1))
    plot.save('cia_'+case+'.pdf')

  