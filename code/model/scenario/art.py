import re
import numpy as np
from copy import deepcopy
from utils import stats,fio,parallel,flatten,minimize,itslice,log,dict_list_update
from model import slicers,params,system,target,fit,out
from model.scenario import tvec,fname,batch_select
import model.scenario

model.scenario.uid = '2022-09-24'
model.scenario.N['sam'] = 100000

cascade = dict( # 2020 cascade targets
  low  = (.60,.40,.80),
  mid  = (.80,.80,.80),
  high = (.95,.95,.95),
)
cases = [
  'fsw-cli+',
  'fsw+cli-',
  'fsw-cli-',
  'fsw+cli+',
]

def get_keyout_data(R,t,case):
  pops = ['all','w','m','aq','wq','mq','fsw','cli']
  to = tvec['outs'].tolist()
  Pkeys = ['seed'] + list(params.def_sample_distrs().keys()) + ['PX_fsw','PX_cli','EHY_acute']
  return dict(case=case,ll=R.get('ll',np.nan),
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
# objective 1: refitting

def run_refit(case,b):
  log(0,'scenario.run_refit: '+case)
  P0s = batch_select(fio.load(fname('npy','fit','Ps',case='base')),b)
  T,PD = get_refit_T_PD(case)
  fun = lambda P: refit_cascade(P,PD,T,tvec['cal'],ftol=.1)
  Ps = parallel.ppool(len(P0s)).map(fun,P0s); log(1)
  fio.save(fname('npy','art','Ps',case=case,b=b),Ps)

def merge_refit():
  log(0,'scenario.merge_refit')
  for case in cases:
    Ps = [P for b in range(model.scenario.N['batch']) for P in
      fio.load(fname('npy','art','Ps',case=case,b=b))]
    fio.save(fname('npy','art','Ps',case=case),Ps)

def rerun_refit():
  log(0,'scenario.rerun_refit')
  for case in ['base']+cases:
    log(1,case)
    T = target.get_all_esw() if case=='base' else get_refit_T_PD(case)[0]
    Ps = fio.load(fname('npy','fit' if case=='base' else 'art','Ps',case=case))
    Rs = system.run_n(Ps,t=tvec['main'],T=T)
    fio.save_csv(fname('csv','art','expo',case=case),
      out.expo(R1s=Rs,tvec=tvec['main'],t=tvec['plot'],snames=['all','w','m','aq','cli','fsw'],
        onames=['incidence','prevalence','diagnosed','treated_c','vls_c','treated_u','vls_u']))
    fio.save_csv(fname('csv','art','keyout',case=case),
      [get_keyout_data(R,tvec['main'],case) for R in Rs])
    fio.save_csv(fname('csv','art','infs',case=case),
      out.get_infections(Rs,tvec['main'],tvec['plot']))
    fit.plot_refit(tvec['main'],Rs,T,fname=fname('fig','art','refit',case=case))
    if case=='base':
      R0s = Rs
      fit.plot_cal(tvec['main'],Rs,T,fname=fname('fig','fit','cal',case=case))
    else:
      fio.save_csv(fname('csv','art','infs-diff',case=case),
        out.get_infections(Rs,tvec['main'],tvec['plot'],R2s=R0s,vsop='1-2'))

def get_refit_T_PD(case):
  parse_case = lambda case: re.findall('(.*?)(\+|\-|\~)',case)
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
  # PD has keys like "rate:pop"; we re-sample based on percentiles (x)
  kwds = dict(method='SLSQP',options=dict(ftol=ftol))
  # we fit relative rates stepwise (dx, +tx, +ux) as it is (usually) more efficient
  PDz = {}; x0z = []; Tz = []; # init step-wise: distrs, x0, targets
  P['refit_iter'] = []
  for pz,oz in zip(['Rdx','Rtx','Rux'],['diagnosed','treated','vls']):
    PDzi = {k:v for k,v in PD.items() if (pz in k)}
    PDz.update(PDzi)
    x0z += [.1]*len(PDzi)
    Tz  += [Ti for Ti in T if (oz in Ti.name)]
    jfun = lambda x: refit_cascade_jfun(x,P,PDz,Tz,t)
    M = minimize(jfun,x0z,bounds=[(0,1) for k in PDz],**kwds)
    x0z[:] = M.x
    if not M.success:
      return False
    else:
      P['refit_iter'].append(M.nfev)
  # DEBUG: non-stepwise version
  # jfun = lambda xs: fit_cascade_jfun(xs,P,PD,T,t)
  # M = minimize(jfun,[.1 for k in PD],bounds=[(0,1) for k in PD],**kwds)
  return P # TODO: is this robust updated P?

def refit_cascade_jfun(x,P,PD,T,t):
  # map xs (quantiles) -> values from PD & update P with the new values
  # e.g. if PD[key] = unif(0,2) and xk = 0.2 -> Pu[key] = 0.4
  P = P_update_Rqx_by_group(P,{key:PD[key].ppf(xk) for key,xk in zip(PD.keys(),x)})
  # run the model & return -ll
  return - system.run(P,t,T,RPts=[])['ll']

def P_update_Rqx_by_group(P,Pu):
  # updadte elements in {Rdx,Rtx,Rux}_scen using Rqx_by_group & kwds in Pu, e.g. Rdx:fsw=1.5
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

def run_sens(N=10):
  P0s = fio.load(fname('npy','fit','Ps',b='all'))
  Ps = get_sens_sample(P0s,N)
  Rs = system.run_n(Ps,t=tvec['main'])
  fio.save_csv(fname('csv','art','keyout',case='sens'),
    [get_keyout_data(R,tvec['main'],'sens') for R in Rs])
  fio.save_csv(fname('csv','art','infs',case='sens'),
    out.get_infections(Rs,tvec['main'],tvec['plot']))
  fit.plot_refit(tvec['main'],Rs,T=None,fname=fname('fig','art','refit',case='sens'))

def get_sens_sample(Ps,N):
  PDs = {
    'd': stats.beta_binom(p=.6,n=3),
    't': stats.beta_binom(p=.6,n=3),
    'u': stats.gamma_p(p=4,v=12),
  }
  Rx_si = np.ones([2,4,1,1])
  PD = {'R'+step+'x:'+pop: PDs[step] for pop in ('fsw','cli','aq') for step in 'dtu' }
  return [P_update_Rqx_by_group(deepcopy(P),Pu)
    for P in dict_list_update(Ps,Rdx_si=Rx_si,Rtx_si=Rx_si)
    for Pu in params.get_n_sample_lhs(N,PD,seed=P['seed'])]

if __name__ == '__main__':
  # rerun_refit()
  run_sens()
