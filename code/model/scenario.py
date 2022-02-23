import re
import numpy as np
from copy import deepcopy
from utils import rootpath,genpath,flatten,minimize,itslice,log
from utils import stats,fio,parallel
from model import slicers,params,system,target,fit,out

# --------------------------------------------------------------------------------------------------
# config stuff

uid = fio.datestamp()
N = dict(
  cal = 100000,
  topcal = .1,
  topfit = .01,
  batch = 10,
  sens = 10,
  b = 0,
)
tvec = dict( # time vectors
  cal = system.f_t(tf=2021),
  fit = system.f_t(tf=2021),
  main = system.f_t(tf=2050,dt=.05),
  infs = system.f_t(tf=2050,dt=1),
  outs = system.f_t(tf=2050,t0=2000,dt=5),
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

# --------------------------------------------------------------------------------------------------
# helpers

def fname(ftype,phase,key,case='base',b=None):
  if b is None: b = N['b']
  subdirs = [uid,str(N['cal'])]
  if ftype=='npy':
    path,ext = ['data','npy',*subdirs],''
  if ftype=='csv':
    path,ext = ['data','mid',*subdirs],'.csv'
  if ftype=='fig':
    path,ext = ['out','fig',*subdirs],'.pdf'
  # e.g. data/npy/2022-01-01/1000/fit_Ps_fsw+cli+_3.npy
  return genpath(rootpath(*path,'{}_{}_{}_{}{}'.format(phase,key,case,b,ext)))

def batch_select(objs):
  nb = int(len(objs) / N['batch'])
  return objs[slice(nb*N['b'],nb*(N['b']+1))]

def parse_case(case):
  # e.g. fsw-cli-aq+  =>  [('fsw','-'),('cli','-'),('aq','+')]
  return re.findall('(.*?)(\+|\-|\~)',case)

def get_expo_data(R,t,case):
  pops = ['all','aq','wq','mq','fsw','cli']
  to = tvec['outs'].tolist()
  Pkeys = ['seed',
    'PX_fsw','PX_cli','dur_fsw_l','dur_fsw_h','dur_cli', # PX & turnover
    'EHY_acute','Rbeta_gud_inf','Rbeta_gud_sus_w','P_gud_fsw_l','iP_gud_cli_fsw:gp', # beta
    'A_mc','A_swr','PA_ai_mcq','PA_ai_swq', # sex acts
    'C_swo_fsw_l','C_swr_fsw_l','RC_swo_fsw_h:l','RC_swr_fsw_h:l','A_swq_cli', # sex work
    'Rdx_fsw','Rdx_cli','Rdx_mqq','Rtx_fsw','Rtx_cli','Rtx_mqq', # diag & treat
  ]
  return dict(case=case,
    **{key:R['P'][key] for key in Pkeys},
    **{'cuminf_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.cuminfect(R,**slicers[pop].pop,tvec=t)[itslice(to,t)])},
    **{'inc_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.incidence(R,**slicers[pop].pop,tvec=t,t=to))},
    **{'prev_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.prevalence(R,**slicers[pop].pop,tvec=t,t=to))},
    **{step+'_'+pop+'_2020':out.by_name(step)(R,**slicers[pop].pop,tvec=t,t=2020)[0] for pop in pops
        for step in ['diagnosed','treated_c','vls_c','treated_u','vls_u']},
  )

# --------------------------------------------------------------------------------------------------
# objective 0: calibration

def run_calibrate():
  log(0,'scenario.run_calibrate')
  T = target.get_all_esw()
  seeds = batch_select(range(N['cal']))
  P0s = params.get_n_all(len(seeds),seeds=seeds)
  R0s = system.run_n(P0s,t=tvec['cal'],T=T)
  R0s = system.drop_fails(R0s)[0]
  fit.plot_cal(tvec['cal'],R0s,T,fname=fname('fig','cal','Ps'))
  Rs = target.top_ll(R0s,top=int(len(seeds)*N['topcal']))
  fio.save(fname('npy','cal','Ps'),[R['P'] for R in Rs])

def merge_calibrate():
  P0s = [P for b in range(N['batch']) for P in fio.load(fname('npy','cal','Ps',b=b))]
  Ps = target.top_ll(P0s,top=int(N['cal']*N['topfit']))
  fio.save(fname('npy','top','Ps',b='all'),Ps)

# --------------------------------------------------------------------------------------------------
# objective 1: refitting

def run_refit():
  log(0,'scenario.run_refit')
  P0s = batch_select(fio.load(fname('npy','top','Ps',b='all')))
  for case in cases:
    log(1,case)
    T,PD = get_refit_T_PD(case)
    fun = lambda P: refit_cascade(P,PD,T,tvec['fit'],ftol=.1)
    Ps = parallel.ppool(len(P0s)).map(fun,P0s); log(1)
    fio.save(fname('npy','fit','Ps',case=case),Ps)

def merge_refit():
  for case in cases:
    Ps = [P for b in range(N['batch']) for P in fio.load(fname('npy','fit','Ps',case=case,b=b))]
    fio.save(fname('npy','fit','Ps',case=case,b='all'),Ps)

def expo_refit():
  log(0,'scenario.expo_refit')
  for case in ['base']+cases:
    log(1,case)
    T = target.get_all_esw() if case=='base' else get_refit_T_PD(case)[0]
    Ps = fio.load(fname('npy','top' if case=='base' else 'fit','Ps',case=case,b='all'))
    Rs = system.run_n(Ps,t=tvec['main'],T=T)
    fio.save_csv(fname('csv','fit','expo',case=case,b='all'),
      [get_expo_data(R,tvec['main'],case) for R in Rs])
    fio.save_csv(fname('csv','fit','infs',case=case,b='all'),
      out.get_infections(Rs,tvec['main'],tvec['infs']))
    fit.plot_refit(tvec['main'],Rs,T,fname=fname('fig','fit','Ps',case=case,b='all'))
    if case=='base':
      R0s = Rs
      fit.plot_cal(tvec['main'],Rs,T,fname=fname('fig','top','Ps',case=case,b='all'))
    else:
      fio.save_csv(fname('csv','fit','infs-diff',case=case,b='all'),
        out.get_infections(Rs,tvec['main'],tvec['infs'],R2s=R0s,vsop='1-2'))

def get_refit_T_PD(case):
  Ts = {
    'fsw-': target.make_targets_2020(cascade['low'], s=0,i=(2,3)),
    'fsw+': target.make_targets_2020(cascade['high'],s=0,i=(2,3),w=1e-6),
    'cli-': target.make_targets_2020(cascade['low'], s=1,i=(2,3)),
    'cli+': target.make_targets_2020(cascade['high'],s=1,i=(2,3),w=1e-6),
    'all-': target.make_targets_2020(cascade['mid'], s=(0,1),i=(0,1,2,3)),
    'all+': target.make_targets_2020(cascade['high'],s=(0,1),i=(0,1,2,3),w=1e-6),
  }
  T = flatten( Ts[pop+c] for pop,c in parse_case(case+'all-') )
  eps = 1e-9
  PDs = {
    'd+': stats.uniform(l=1-eps,h=1+eps), 'd-': stats.uniform(l=0,h=1),
    't+': stats.uniform(l=1-eps,h=1+eps), 't-': stats.uniform(l=0,h=1),
    'u+': stats.uniform(l=1-eps,h=1+eps), 'u-': stats.uniform(l=1,h=25),
  }
  PD = {'R'+step+'x:'+pop: PDs[step+c] for pop,c in parse_case(case+'aq-') for step in 'dtu'}
  return T,PD

def refit_cascade(P,PD,T,t,ftol=.1):
  # fit the model to "T", with base parameters P, and re-sampled values in PD
  # PD has keys like "rate:pop"; we re-sample based on percentiles (pfit)
  kwds = dict(method='SLSQP',options=dict(ftol=ftol))
  # we fit relative rates stepwise (dx, +tx, +ux) as it is (usually) more efficient
  PDz = {}; P0z = []; Tz = []; # init step-wise: distrs, x0, targets
  P['refit_iter'] = []
  for pz,oz in zip(['Rdx','Rtx','Rux'],['diagnosed','treated','vls']):
    PDzi = {k:v for k,v in PD.items() if (pz in k)}
    PDz.update(PDzi)
    P0z += [.1 for k in PDzi]
    Tz  += [Ti for Ti in T if (oz in Ti.name)]
    jfun = lambda pfits: refit_cascade_jfun(pfits,P,PDz,Tz,t)
    M = minimize(jfun,P0z,bounds=[(0,1) for k in PDz],**kwds)
    P0z[:] = M.x
    if not M.success:
      return False
    else:
      P['refit_iter'].append(M.nfev)
  # DEBUG: non-stepwise version
  # jfun = lambda pfits: fit_cascade_jfun(pfits,P,PD,T,t)
  # M = minimize(jfun,[.1 for k in PD],bounds=[(0,1) for k in PD],**kwds)
  return P # TODO: is this robust updated P?

def refit_cascade_jfun(pfits,P,PD,T,t):
  # map pfits (quantiles) -> values from PD & update P with the new values
  P = P_update_Rqx_by_group(P,{key:PD[key].ppf(pfit) for key,pfit in zip(PD.keys(),pfits)})
  # run the model & return -ll
  return - system.run(P,t,T,RPts=[])['ll']

def P_update_Rqx_by_group(P,Pu):
  for key in Pu.keys():
    rate,pop = key.split(':')
    P.update({rate+'_scen':Rqx_by_group(Rqx=P.get(rate+'_scen'),**{pop:Pu[key]})})
  return P

def Rqx_by_group(Rqx=None,**kwds):
  # create/update a relative-rate using pop=value kwds, e.g. FSW=1.5
  if Rqx is None: Rqx = np.ones((2,4,1))
  for sname,sRqx in kwds.items():
    for sk in flatten(slicers[sname].pop.get('s',(0,1))):
      for ik in flatten(slicers[sname].pop.get('i',(0,1,2,3))):
        Rqx[sk,ik,0] = sRqx
  return Rqx

# --------------------------------------------------------------------------------------------------
# objective 2: sensitivity

def run_sens(expo=True,infs=True,plot=True):
  P0s = batch_select(fio.load(fname('npy','top','Ps',b='all')))
  Ps = get_sens_sample(P0s)
  fio.save(fname('npy','sens','Ps'),Ps)
  Rs = system.run_n(Ps,t=tvec['main'])
  # TODO: clean up
  if infs: fio.save_csv(fname('csv','sens','infs',case='sens'),out.get_infections(Rs,tvec['main'],tvec['infs']))
  if expo: fio.save_csv(fname('csv','sens','expo',case='sens'),[get_expo_data(R,tvec['main'],'sens') for R in Rs])
  if plot: fit.plot_refit(tvec['main'],Rs,T=None,fname=fname('fig','sens','Ps',case='sens'))

def get_sens_sample(Ps):
  PDs = {
    'd': stats.beta_binom(p=.6,n=3),
    't': stats.beta_binom(p=.6,n=3),
    'u': stats.gamma_p(p=4,v=12),
  }
  PD = {'R'+step+'x:'+pop: PDs[step] for pop in ('fsw','cli','aq') for step in 'dtu' }
  return [P_update_Rqx_by_group(deepcopy(P),Pu)
    for P in Ps for Pu in params.get_n_sample_lhs(N['sens'],PD,seed=P['seed'])]
