import re,time
import numpy as np
from copy import deepcopy
from utils import rootpath,genpath,flatten,minimize,itslice,squarish
from utils import stats,fio,deco,parallel
from model import slicers,params,system,target,out,plot,fit

plotsize = 3 # inches
uid = fio.datestamp()
cascade_low  = (.40,.60,.80)
cascade_mid  = (.60,.80,.80)
cascade_high = (.95,.95,.95)
cases = [
  'FSW-Cli+',
  'FSW+Cli-',
  'FSW-Cli-',
  'FSW+Cli+',
]

def fname(ftype,key,case,N,N0):
  # TODO: omit N0 if N0==0
  if ftype=='npy':
    path,ext = ['data','.npy',uid],''
  if ftype=='csv':
    path,ext = ['data','mid',uid],'.csv'
  if ftype=='fig':
    path,ext = ['out','fig',uid],'.pdf'
  return genpath(rootpath(*path,'{}_{}_N={}-{}{}'.format(key,case,N0,N0+N-1,ext)))

def sample_run_base(t,T,N,N0=0,top=.10):
  Ps_sam = params.get_n_all(N,seeds=range(N0,N0+N))
  Rs_sam = system.run_n(Ps_sam,t,T)
  Rs_sam = system.drop_fails(Rs_sam)[0]
  fit.plot_all(t,Rs_sam,T,fname=fname('fig','fit','sam',N,N0))
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

def sample_random_lower(Ps,n=1):
  PD = get_refit_case('FSW-Cli-')[1]
  Pss = [P_update_Rqx_by_group(deepcopy(P),Pu)
    for P in Ps for Pu in params.get_n_sample_lhs(n,PD,seed=P['seed'])]
  return Pss

def get_sens_data(R,t,case):
  # extract outputs of interest & selected model params for sensitivity analysis
  to  = list(range(2000,2050+1,5))
  ito = itslice(to,t)
  Pkeys = ['seed',
    'PX_fsw','PX_cli','dur_fsw_l','dur_fsw_h','dur_cli', # PX & turnover
    'EHY_acute','Rbeta_gud_inf','Rbeta_gud_sus_w','P_gud_fsw_l','iP_gud_cli_fsw:gp', # beta
    'A_mc','A_reg','PA_ai_mcq','PA_ai_swq', # sex acts
    'C_new_fsw_l','C_reg_fsw_l','RC_new_fsw_h:l','RC_reg_fsw_h:l','A_swq_cli', # sex work
    'Rdx_fsw','Rdx_cli','Rdx_mqq','Rtx_fsw','Rtx_cli','Rtx_mqq', # diag & treat
  ]
  return dict(case=case,
    **{key:R['P'][key] for key in Pkeys},
    # TODO: add group-specific HIV prevalence
    **{'cuminf_'+str(toi):oi for toi,oi in zip(to,out.cuminfect(R,**slicers['ALL'].pop,tvec=t)[ito])},
    **{'incid_' +str(toi):oi for toi,oi in zip(to,out.incidence(R,**slicers['ALL'].pop,tvec=t,t=to))},
    **{step+'_'+pop:out.by_name(step)(R,**slicers[pop].pop,tvec=t,t=2020)[0]
        for step in ['diagnosed','treated_c','vls_c','treated_u','vls_u']
        for pop in ['ALL','AQ','FSW','Cli']},
    **{
      'inc_all':  out.incidence(R,**slicers['ALL'].pop,tvec=t,t=2020)[0],
      'prev_all': out.prevalence(R,**slicers['ALL'].pop,tvec=t,t=2020)[0],
      'prev_fsw': out.prevalence(R,**slicers['FSW'].pop,tvec=t,t=2020)[0],
      'prev_cli': out.prevalence(R,**slicers['Cli'].pop,tvec=t,t=2020)[0],
      'prev_aq':  out.prevalence(R,**slicers['AQ'].pop,tvec=t,t=2020)[0],
      'prev_ratio_fsw.wq': out.vs_pop('prevalence',R,slicers['FSW'].pop,slicers['WQ'].pop,tvec=t,t=2020)[0],
      'prev_ratio_cli.mq': out.vs_pop('prevalence',R,slicers['Cli'].pop,slicers['MQ'].pop,tvec=t,t=2020)[0],
    }
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
    **summary('prev_all',[out.prevalence(R,**slicers['ALL'].pop) for R in Rs]),
    **summary('inc_all', [out.incidence(R,**slicers['ALL'].pop) for R in Rs]),
    **summary('inc_aq',  [out.incidence(R,**slicers['AQ'].pop)  for R in Rs]),
    **summary('inc_fsw', [out.incidence(R,**slicers['FSW'].pop) for R in Rs]),
    **summary('inc_cli', [out.incidence(R,**slicers['Cli'].pop) for R in Rs]),
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
  # e.g. FSW-Cli-AQ+  =>  [('FSW','-'),('Cli','-'),('AQ','+')]
  return re.findall('(.*?)(\+|\-)',case)

def get_refit_case(case):
  eps = 1e-9
  PDs = {
    'd+': stats.uniform(l=1-eps,h=1+eps), 'd-': stats.uniform(l=0,h=1),
    't+': stats.uniform(l=1-eps,h=1+eps), 't-': stats.uniform(l=0,h=1),
    'u+': stats.uniform(l=1-eps,h=1+eps), 'u-': stats.uniform(l=1,h=25),
  }
  T2s = {
    'FSW-': target.make_targets_2020(cascade_low, s=0,i=(2,3)),
    'FSW+': target.make_targets_2020(cascade_high,s=0,i=(2,3),w=1e-6),
    'Cli-': target.make_targets_2020(cascade_low, s=1,i=(2,3)),
    'Cli+': target.make_targets_2020(cascade_high,s=1,i=(2,3),w=1e-6),
    'ALL-': target.make_targets_2020(cascade_mid, s=(0,1),i=(0,1,2,3)),
    'ALL+': target.make_targets_2020(cascade_high,s=(0,1),i=(0,1,2,3),w=1e-6),
  }
  T2 = flatten( T2s[pop+c] for pop,c in parse_case(case+'ALL-') )
  PD = {'R'+step+'x:'+pop: PDs[step+c] for pop,c in parse_case(case+'AQ-') for step in 'dtu'}
  return T2,PD

def do_scenario(Ps,t2,T2,case,N,N0,outs=True,sens=True,infs=True,plotfit=True):
  Rs = system.run_n(Ps,t2,T2)
  if outs: fio.save_csv(fname('csv','outs',case,N,N0),get_outs(Rs,t2,case))
  if sens: fio.save_csv(fname('csv','sens',case,N,N0),get_sens_data_n(Rs,t2,case))
  if infs: fio.save_csv(fname('csv','infs',case,N,N0),out.get_infections(Rs,t2,system.f_t(dt=1)))
  if plotfit: fit.plot_refit(t2,Rs,T2,fname('fig','fit',case,N,N0))
  return Rs

def main_fit(N,N0=0,sample=True,refit=True,top=.10,outs=True,sens=True,infs=True,plotfit=True):
  t  = system.f_t(t1=2021)
  t2 = system.f_t(t1=2050,dt=.05)
  T1 = target.get_all_esw()
  case = 'base'
  if sample:
    R0s = sample_run_base(t,T1,N,N0=N0,top=top)
    P0s = fio.save(fname('npy','Ps',case,N,N0),[R['P'] for R in R0s])
  else:
    P0s = fio.load(fname('npy','Ps',case,N,N0))
  R0s = do_scenario(P0s,t2,T1,case,N,N0,outs=outs,sens=sens,infs=infs,plotfit=plotfit)
  for case in cases: # counterfactuals (fitted)
    print('\n'+case,flush=True)
    T2,PD = get_refit_case(case)
    if refit:
      P1s = fit_cascade_n(P0s,PD,T2,t,ftol=.1)
      fio.save(fname('npy','Ps',case,N,N0),P1s)
    else:
      P1s = fio.load(fname('npy','Ps',case,N,N0))
    # re-run for t2
    R1s = do_scenario(P1s,t2,T2,case,N,N0,outs=outs,sens=sens,infs=infs,plotfit=plotfit)
    if infs: fio.save_csv(fname('csv','infs-diff',case,N,N0),
      out.get_infections(R1s,t2,system.f_t(dt=1),R2s=R0s,vsop='1-2'))
    # plot cumulative additional infections
    plot_diff(t2,R1s,R0s,'cuminfect',['ALL','AQ','FSW','Cli'],intervals=[.5,.95],ylim=(0,1),
      xlim=(1998,2042),vsop='1-2/2',ylab='Cumulative Additional Infections')
    plot.save(fname('fig','cai4',case,N,N0))
    if N > 256: time.sleep(30) # let CPU cool :)

def main_random(N,N0=0,Nr=10,sample=True,resample=True,top=.10,outs=True,sens=True,infs=True,plotfit=True):
  t  = system.f_t(t1=2021)
  t2 = system.f_t(t1=2050,dt=.05)
  T1 = target.get_all_esw()
  # base case (fitted)
  case = 'base'
  if sample:
    R0s = sample_run_base(t,T1,N,N0=N0,top=top)
    P0s = fio.save(fname('npy','Ps',case,N,N0),[R['P'] for R in R0s])
  else:
    P0s = fio.load(fname('npy','Ps',case,N,N0))
  R0s = do_scenario(P0s,t2,T1,case,N,N0,outs=outs,sens=sens,infs=infs,plotfit=plotfit)
  # counterfactuals (random)
  case = 'RL'+str(Nr)
  if resample:
    P1s = sample_random_lower(P0s,Nr)
    fio.save(fname('npy','Ps',case,N,N0),P1s)
  else:
    P1s = fio.load(fname('npy','Ps',case,N,N0))
  R1s = do_scenario(P1s,t2,None,case,N,N0,outs=outs,sens=sens,infs=infs,plotfit=plotfit)
