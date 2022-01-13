import numpy as np
from utils import stats,rootpath,genpath,flatten,ppool,minimize,itslice
from utils import fio,deco
from model import slicers,params,system,target,out,plot,handfit

plotsize = 3 # inches
uid = fio.datestamp()

def fname(case,key,N,N0):
  if case=='npy':
    path,ext = ['data','.npy',uid],''
  if case=='csv':
    path,ext = ['data','mid',uid],'.csv'
  if case=='fig':
    path,ext = ['out','fig',uid],'.pdf'
  return genpath(rootpath(*path,'{}_N={}-{}{}'.format(key,N0,N0+N-1,ext)))

def get_sample(t,T,N,N0=0,top=.10):
  Ps_sam = params.get_n_all(N,seeds=range(N0,N0+N))
  Rs_sam = system.run_n(Ps_sam,t,T)
  Rs_sam = system.drop_fails(Rs_sam)[0]
  handfit.plot_all(t,Rs_sam,T,fname=fname('fig','handfit_sam',N,N0))
  Rs = target.top_q_ll(Rs_sam,top=top)
  return Rs

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

def get_sens_data(R,t,case):
  # extract outputs of interest & selected model params for sensitivity analysis
  to  = [2020,2025,2030,2035,2040,2045,2050]
  ito = itslice(to,t)
  Pkeys = ['seed',
    'beta_0','EHM_acute','Rbeta_vr', # HIV
    'A_mc','RPA_condom_av','PA_circum_2050','PA_ai_mc','PA_ai_sw',# sex
    'RP_gud_2050','Rbeta_gud_s','Rbeta_gud_i' # GUD
    'PX_fsw','dur_fsw_l','C_new_fswl', # FSW
    'A_swq_cli','dur_cli', # Cli
  ]
  D = dict(case=case,
    **{'cuminf_'+str(toi):oi for toi,oi in zip(to,out.cuminfect(R,**slicers['ALL'].pop,aggr=True,tvec=t)[ito])},
    **{'incid_' +str(toi):oi for toi,oi in zip(to,out.incidence(R,**slicers['ALL'].pop,aggr=True,tvec=t,t=to))},
    **{key:R['P'][key] for key in Pkeys},
  )
  return D

def get_sens_data_n(Rs,t,case):
  return [get_sens_data(R,t,case) for R in Rs]

def plot_diff(t,R1s,R2s,oname,snames,intervals=.95,ylim=None,**kwds):
  assert len(R1s)==len(R2s), 'R1s and R2s must have the same length'
  R1s,R2s = system.drop_fails(R1s,R2s)
  fh,ah = plot.subplots(1,len(snames))
  intervals = flatten(intervals)
  for i,sname in enumerate(snames): # subplots
    plot.plt.sca(ah[0,i])
    plot.labels(title=slicers[sname].label,
      x='Time (years)',y=out.labels.get(oname,oname) if i==0 else None)
    for interval in intervals:
      plot.plot_SvR(oname,t,R1s,R2s,sname,vsop='1-2/1',tvec=t,interval=interval)
    plot.plt.ylim(ylim)
  fh.set_size_inches((plotsize*len(snames),plotsize))
  fh.tight_layout()
  return fh,ah

def get_refit_case(case):
  eps = 1e-9
  PDD = {
    '0-1':  stats.uniform(l=0,h=1),
    '1-1':  stats.uniform(l=1-eps,h=1+eps),
    '1-50': stats.uniform(l=1,h=50),
  }
  if case=='LoLo':
    T2 = target.get_cascade_2020('ssa-lowest')\
       + target.get_cascade_2020('ssa-lowest-fsw')
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
  return T2,PD

def main(N,N0=0,sample=True,refit=True,top=.10,sens=True):
  t  = system.f_t(t1=2025)
  t2 = system.f_t(t1=2050.05,dt=.05) # TEMP
  T1 = target.get_all_esw()
  # base case (fitted)
  if sample:
    R1s = get_sample(t,T1,N,N0=N0,top=top)
    P1s = fio.save(fname('npy','P1s',N,N0),[R['P'] for R in R1s])
  else:
    P1s = fio.load(fname('npy','P1s',N,N0))
  # re-run for t2
  R1s = system.run_n(P1s,t2,T1)
  if sens: fio.save_csv(fname('csv','sens_base',N,N0),get_sens_data_n(R1s,t2,'base'))
  handfit.plot_all(t2,R1s,T1,fname=fname('fig','handfit_base',N,N0)) # DEBUG
  # counterfactuals
  for case in ['LoLo','LoHi','HiLo']:
    T2,PD = get_refit_case(case)
    if refit:
      P2s = fit_cascade_n(P1s,PD,T2,t,ftol=.01)
      fio.save(fname('npy','P2s_'+case,N,N0),P2s)
    else:
      P2s = fio.load(fname('npy','P2s_'+case,N,N0))
    # re-run for t2
    R2s = system.run_n(P2s,t2,T2)
    if sens: fio.save_csv(fname('csv','sens_'+case,N,N0),get_sens_data_n(R2s,t2,case))
    handfit.plot_all(t2,R2s,T2,fname=fname('fig','handfit_'+case,N,N0)) # DEBUG
    # plot cia: cumulative infections averted
    plot_diff(t2,R2s,R1s,'cuminfect',['ALL','FSW','Cli'],intervals=[.5,.95],ylim=(0,1))
    plot.save(fname('fig','cia_'+case,N,N0))
