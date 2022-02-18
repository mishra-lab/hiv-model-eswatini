import re,time
import numpy as np
from copy import deepcopy
from utils import rootpath,genpath,flatten,minimize,itslice,squarish
from utils import stats,fio,deco,parallel
from model import slicers,params,system,target,out,plot,fit
plotsize = 3 # inches

# override any of these global config values in main before running
uid = fio.datestamp() # unique ID
N = dict(
  top = .10,
  size = 1000,
  batch = 0,
  rand = 10,
)
tvec = dict( # time vectors
  fits = system.f_t(tf=2021),
  main = system.f_t(tf=2050,dt=.05),
  infs = system.f_t(tf=2050,dt=1),
) 
cascade = dict( # 2020 cascade targets
  low  = (.40,.60,.80),
  mid  = (.60,.80,.80),
  high = (.95,.95,.95),
)
cases = [
  'fsw-cli+',
  'fsw+cli-',
  'fsw-cli-',
  'fsw+cli+',
]

def fname(ftype,key,case):
  if ftype=='npy':
    path,ext = ['data','.npy',uid],''
  if ftype=='csv':
    path,ext = ['data','mid',uid],'.csv'
  if ftype=='fig':
    path,ext = ['out','fig',uid],'.pdf'
  return genpath(rootpath(*path,'{}_{}_N={}-{}{}'.format(key,case,N['size'],N['batch'],ext)))

def sample_run_base(t,T):
  seeds = range(N['batch']*N['size'],(N['batch']+1)*N['size'])
  Ps_sam = params.get_n_all(N['size'],seeds=seeds)
  Rs_sam = system.run_n(Ps_sam,t,T)
  Rs_sam = system.drop_fails(Rs_sam)[0]
  fit.plot_all(t,Rs_sam,T,fname=fname('fig','fit','sam'))
  Rs = target.top_q_ll(Rs_sam,top=N['top'])
  return Rs

def Rqx_by_group(Rqx=None,**kwds):
  # create/update a relative-rate using pop=value kwds, e.g. FSW=1.5
  if Rqx is None: Rqx = np.ones((2,4,1))
  for sname,sRqx in kwds.items():
    for sk in flatten(slicers[sname].pop.get('s',(0,1))):
      for ik in flatten(slicers[sname].pop.get('i',(0,1,2,3))):
        Rqx[sk,ik,0] = sRqx
  return Rqx

def P_update_Rqx_by_group(P,Pu):
  for key in Pu.keys():
    rate,pop = key.split(':')
    P.update({rate+'_scen':Rqx_by_group(Rqx=P.get(rate+'_scen'),**{pop:Pu[key]})})
  return P

def fit_cascade_jfun(pfits,P,PD,T,t):
  # map pfits (quantiles) -> values from PD & update P with the new values
  P = P_update_Rqx_by_group(P,{key:PD[key].ppf(pfit) for key,pfit in zip(PD.keys(),pfits)})
  # run the model
  return - system.run(P,t,T,RPts=[])['ll']

def fit_cascade(P,PD,T,t,ftol=.1):
  # TODO: re-implement specified P0
  # fit the model to "T", with base parameters P, and re-sampled values in PD
  # PD has keys like "rate:pop"; we re-sample based on percentiles (pfit)
  kwds = dict(method='SLSQP',options=dict(ftol=ftol))
  # we fit relative rates stepwise (dx, +tx, +ux) as it is more efficient
  PDz = {}; P0z = []; Tz = []; # init step-wise: distrs, x0, targets
  P['refit_iter'] = []
  for pz,oz in zip(['Rdx','Rtx','Rux'],['diagnosed','treated','vls']):
    PDzi = {k:v for k,v in PD.items() if (pz in k)}
    PDz.update(PDzi)
    P0z += [.1 for k in PDzi]
    Tz  += [Ti for Ti in T if (oz in Ti.name)]
    jfun = lambda pfits: fit_cascade_jfun(pfits,P,PDz,Tz,t)
    M = minimize(jfun,P0z,bounds=[(0,1) for k in PDz],**kwds)
    P0z[:] = M.x
    if not M.success:
      return False
    else:
      P['refit_iter'].append(M.nfev)
  # DEBUG: non-stepwise
  # jfun = lambda pfits: fit_cascade_jfun(pfits,P,PD,T,t)
  # M = minimize(jfun,[.1 for k in PD],bounds=[(0,1) for k in PD],**kwds)
  return P # TODO: is this robust updated P?

def fit_cascade_n(Ps,PD,T,t,**kwds):
  fun = lambda P: fit_cascade(P,PD=PD,T=T,t=t,**kwds)
  return parallel.ppool(len(Ps)).map(fun,Ps)

def sample_random_lower(Ps):
  PD = get_refit_case('fsw~cli~aq~',default=False)[1]
  Pss = [P_update_Rqx_by_group(deepcopy(P),Pu)
    for P in Ps for Pu in params.get_n_sample_lhs(N['rand'],PD,seed=P['seed'])]
  return Pss

def get_sens_data(R,t,case):
  # extract outputs of interest & selected model params for sensitivity analysis
  pops = ['all','aq','fsw','cli']
  to  = list(range(2000,2050+1,5))
  ito = itslice(to,t)
  Pkeys = ['seed',
    'PX_fsw','PX_cli','dur_fsw_l','dur_fsw_h','dur_cli', # PX & turnover
    'EHY_acute','Rbeta_gud_inf','Rbeta_gud_sus_w','P_gud_fsw_l','iP_gud_cli_fsw:gp', # beta
    'A_mc','A_swr','PA_ai_mcq','PA_ai_swq', # sex acts
    'C_swo_fsw_l','C_swr_fsw_l','RC_swo_fsw_h:l','RC_swr_fsw_h:l','A_swq_cli', # sex work
    'Rdx_fsw','Rdx_cli','Rdx_mqq','Rtx_fsw','Rtx_cli','Rtx_mqq', # diag & treat
  ]
  return dict(case=case,
    **{key:R['P'][key] for key in Pkeys},
    # TODO: add group-specific HIV prevalence
    **{'cuminf_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.cuminfect(R,**slicers[pop].pop,tvec=t)[ito])},
    **{'inc_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.incidence(R,**slicers[pop].pop,tvec=t,t=to))},
    **{'prev_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.prevalence(R,**slicers[pop].pop,tvec=t,t=to))},
    **{'prev_ratio_'+pop[0]+'.'+pop[1]+'_'+str(toi):oi for pop in [('fsw','wq'),('cli','mq')]
        for toi,oi in zip(to,out.vs_pop('prevalence',R,slicers[pop[0]].pop,slicers[pop[1]].pop,tvec=t,t=to))},
    **{step+'_'+pop+'_2020':out.by_name(step)(R,**slicers[pop].pop,tvec=t,t=2020)[0] for pop in pops
        for step in ['diagnosed','treated_c','vls_c','treated_u','vls_u']},
  )

def get_sens_data_n(Rs,t,case):
  return [get_sens_data(R,t,case) for R in Rs]

def get_outs(Rs,t,case):
  summary = lambda oname,Os: {
    oname+'.mu': np.mean(Os,axis=0),
    oname+'.md': np.median(Os,axis=0),
    oname+'.lo': np.quantile(Os,.025,axis=0),
    oname+'.q1': np.quantile(Os,.25, axis=0),
    oname+'.q3': np.quantile(Os,.75, axis=0),
    oname+'.hi': np.quantile(Os,.975,axis=0),
  }
  return dict(case=[case]*len(t),t=t,
    **summary('prev_all',[out.prevalence(R,**slicers['all'].pop) for R in Rs]),
    **summary('inc_all', [out.incidence(R,**slicers['all'].pop) for R in Rs]),
    **summary('inc_aq',  [out.incidence(R,**slicers['aq'].pop)  for R in Rs]),
    **summary('inc_fsw', [out.incidence(R,**slicers['fsw'].pop) for R in Rs]),
    **summary('inc_cli', [out.incidence(R,**slicers['cli'].pop) for R in Rs]),
  )

def plot_diff(t,R1s,R2s,oname,snames,vsop,intervals=.95,ylim=None,xlim=None,ylab=None,**kwds):
  assert len(R1s)==len(R2s), 'R1s and R2s must have the same length'
  R1s,R2s = system.drop_fails(R1s,R2s)
  fh,ah = plot.subplots(1,len(snames))
  if ylab is None: ylab = out.labels.get(oname,oname)
  for i,sname in enumerate(snames): # subplots
    plot.plt.sca(ah[0,i])
    plot.labels(title=slicers[sname].label,x='Time (years)',y=ylab if i==0 else None)
    for interval in flatten(intervals):
      plot.plot_SvR(oname,t,R1s,R2s,sname,vsop=vsop,tvec=t,box=5,**kwds)
    plot.lims(y=ylim,x=xlim)
    plot.plt.grid(axis='y',lw=.5,zorder=0,color=(.8,.8,.8))
  fh.set_size_inches((plotsize*len(snames),plotsize))
  fh.tight_layout()
  return fh,ah

def parse_case(case):
  # e.g. fsw-cli-aq+  =>  [('fsw','-'),('cli','-'),('aq','+')]
  return re.findall('(.*?)(\+|\-|\~)',case)

def get_refit_case(case,default=True):
  eps = 1e-9
  case_PD = case+'aq-'  if default else case
  case_T2 = case+'all-' if default else case
  PDs = {
    'd+': stats.uniform(l=1-eps,h=1+eps), 'd-': stats.uniform(l=0,h=1),  'd~': stats.beta_binom(p=.6,n=3),
    't+': stats.uniform(l=1-eps,h=1+eps), 't-': stats.uniform(l=0,h=1),  't~': stats.beta_binom(p=.6,n=3),
    'u+': stats.uniform(l=1-eps,h=1+eps), 'u-': stats.uniform(l=1,h=25), 'u~': stats.gamma_p(p=4,v=12),
  }
  T2s = {
    'fsw-': target.make_targets_2020(cascade['low'], s=0,i=(2,3)),
    'fsw+': target.make_targets_2020(cascade['high'],s=0,i=(2,3),w=1e-6),
    'cli-': target.make_targets_2020(cascade['low'], s=1,i=(2,3)),
    'cli+': target.make_targets_2020(cascade['high'],s=1,i=(2,3),w=1e-6),
    'all-': target.make_targets_2020(cascade['mid'], s=(0,1),i=(0,1,2,3)),
    'all+': target.make_targets_2020(cascade['high'],s=(0,1),i=(0,1,2,3),w=1e-6),
    'fsw~': [], 'cli~': [], 'aq~': [], # TODO: clean this up?
  }
  T2 = flatten( T2s[pop+c] for pop,c in parse_case(case_T2) )
  PD = {'R'+step+'x:'+pop: PDs[step+c] for pop,c in parse_case(case_PD) for step in 'dtu'}
  return T2,PD

def do_scenario(Ps,t,T2,case,outs=True,sens=True,infs=True,plotfit=True):
  Rs = system.run_n(Ps,t,T2)
  if outs: fio.save_csv(fname('csv','outs',case),get_outs(Rs,t,case))
  if sens: fio.save_csv(fname('csv','sens',case),get_sens_data_n(Rs,t,case))
  if infs: fio.save_csv(fname('csv','infs',case),out.get_infections(Rs,t,tvec['infs']))
  if plotfit: fit.plot_refit(t,Rs,T2,fname('fig','fit',case))
  return Rs

def main_fit(sample=True,refit=True,outs=True,sens=True,infs=True,plotfit=True):
  T1 = target.get_all_esw()
  case = 'base'
  if sample:
    R0s = sample_run_base(tvec['fits'],T1)
    P0s = fio.save(fname('npy','Ps',case),[R['P'] for R in R0s])
  else:
    P0s = fio.load(fname('npy','Ps',case))
  R0s = do_scenario(P0s,tvec['main'],T1,case,outs=outs,sens=sens,infs=infs,plotfit=plotfit)
  for case in cases: # counterfactuals (fitted)
    print('\n'+case,flush=True)
    T2,PD = get_refit_case(case)
    if refit:
      P1s = fit_cascade_n(P0s,PD,T2,tvec['fits'],ftol=.1)
      fio.save(fname('npy','Ps',case),P1s)
    else:
      P1s = fio.load(fname('npy','Ps',case))
    # re-run for t2
    R1s = do_scenario(P1s,tvec['main'],T2,case,outs=outs,sens=sens,infs=infs,plotfit=plotfit)
    if infs: fio.save_csv(fname('csv','infs-diff',case),
      out.get_infections(R1s,tvec['main'],tvec['infs'],R2s=R0s,vsop='1-2'))
    # plot cumulative additional infections
    plot_diff(tvec['main'],R1s,R0s,'cuminfect',['all','aq','fsw','cli'],intervals=[.5,.95],ylim=(0,1),
      xlim=(1998,2042),vsop='1-2/2',ylab='Cumulative Additional Infections')
    plot.save(fname('fig','cai4',case))
    if refit and N['size'] > 256: time.sleep(30) # let CPU cool :)

def main_random(sample=True,resample=True,outs=True,sens=True,infs=True,plotfit=True):
  T1 = target.get_all_esw()
  # base case (fitted)
  case = 'base'
  if sample:
    R0s = sample_run_base(tvec['fits'],T1)
    P0s = fio.save(fname('npy','Ps',case),[R['P'] for R in R0s])
  else:
    P0s = fio.load(fname('npy','Ps',case))
  R0s = do_scenario(P0s,tvec['main'],T1,case,outs=outs,sens=sens,infs=infs,plotfit=plotfit)
  # counterfactuals (random)
  case = 'RL'+str(N['rand'])
  if resample:
    P1s = sample_random_lower(P0s)
    fio.save(fname('npy','Ps',case),P1s)
  else:
    P1s = fio.load(fname('npy','Ps',case))
  R1s = do_scenario(P1s,tvec['main'],None,case,outs=outs,sens=sens,infs=infs,plotfit=plotfit)
