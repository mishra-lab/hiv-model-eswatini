import numpy as np
from copy import deepcopy
from utils import stats,log,fio,ppool,flatten,minimize
from model import strats,params,system,target,fit,out
from model.scenario import tvec,fname

cascade = dict( # 2020 cascade targets
  low  = (.60,.40,.80),
  mid  = (.80,.80,.90),
  high = (.95,.95,.95),
)
cases = ['fsw+cli+','fsw+cli-','fsw-cli+','fsw-cli-']

tkp = dict(tvec=tvec['main'],t=tvec['plot'])
tko = dict(tvec=tvec['main'],t=tvec['outs'])
ekwds = dict(
  skeys = ['all','w','m','aq','fsw','cli'],
  onames = ['incidence','prevalence','cuminfect','diagnosed','treated_c','treated_u','vls_c','vls_u'],
  mode = 'id')

def run_rf(b):
  case = cases[b]
  log(0,'art.run_rf: {}'.format(case))
  P0s = fio.load_npy(fname('npy','fit','Ps',case='base'))
  T = get_refit_T(case+'all-')
  fun = lambda P: run_rf_1(P,T,tvec['cal'])
  Ps = ppool().map(fun,P0s); log(1)
  fio.save_npy(fname('npy','art-rf','Ps',case=case),Ps)

def rerun_rf():
  log(0,'art.rerun_rf')
  for case in ['base']+cases:
    log(1,case)
    base = (case == 'base')
    T = get_refit_T('all+' if base else case+'all-')
    Ps = fio.load_npy(fname('npy','fit' if base else 'art-rf','Ps',case=case))
    R1s = system.run_n(Ps,t=tvec['main'],T=T)
    fio.save_csv(fname('csv','art-rf','wiw',case=case),out.wiw(R1s,**tkp))
    fio.save_csv(fname('csv','art-rf','expo',case=case),out.expo(R1s,**tkp,**ekwds))
    fit.plot_sets(tvec['main'],R1s,T=T,tfname=fname('fig','art-rf','{}',case=case),
      sets='cascade',skeys=['all','aq','fsw','cli'])
    if base:
      fio.save_csv(fname('csv','art-ss','expo',case=case),out.expo(R1s,**tko,**ekwds))
      R2s = deepcopy(R1s)
    else:
      fio.save_csv(fname('csv','art-rf','wiw-diff',case=case),out.wiw(R1s,R2s=R2s,vsop='1-2',**tkp))

def get_refit_T(case):
  T = {
    'fsw-': target.make_targets_2020(cascade['low'], s=0,i=(2,3)),
    'fsw+': target.make_targets_2020(cascade['mid'], s=0,i=(2,3)),
    'cli-': target.make_targets_2020(cascade['low'], s=1,i=(2,3)),
    'cli+': target.make_targets_2020(cascade['mid'], s=1,i=(2,3)),
    'all-': target.make_targets_2020(cascade['mid'], s=(0,1),i=(0,1,2,3)),
    'all+': target.make_targets_2020(cascade['high'],s=(0,1),i=(0,1,2,3)),
  }
  return flatten( T[k] for k in T if k in case )

def run_rf_1(P,T,t):
  Z0 = [
    dict(key='Rdx',D=stats.uniform(l=0,h= 1),x0=[.5,.5,.5],out='diagnosed'),
    dict(key='Rtx',D=stats.uniform(l=0,h= 1),x0=[.5,.5,.5],out='treated'),
    dict(key='Rux',D=stats.uniform(l=1,h=20),x0=[.5,.5,.5],out='vls'),
  ]
  def run_rf_step(Z,tol,x0=None,**kwds):
    log(3,'[{}]'.format('.'.join([Zi['key'][1] for Zi in Z])).rjust(9)+' ')
    Dz = {Zi['key']+':'+skey: Zi['D'] for Zi in Z for skey in ('fsw','cli','aq')}
    Tz = flatten(Ti for Zi in Z for Ti in T if (Zi['out'] in Ti.name))
    xz = flatten(Zi['x0'] for Zi in Z) if x0 is None else x0
    jfun = lambda x: - system.run(Rxs_update(P,Dz,q=x,**kwds),t,T,RPts=[])['ll']
    M = minimize(jfun,xz,bounds=[(0,1) for x in xz],method='L-BFGS-B',options=dict(ftol=tol))
    Rxs_update(P,Dz,q=M.x,**kwds)
    return M.x, {key:P[key] for Zi in Z if (key := Zi['key']+'_scen')}
  xd,kwd = run_rf_step([Z0[0]],tol=1)
  xt,kwt = run_rf_step([Z0[1]],tol=1,**kwd)
  xu,kwu = run_rf_step([Z0[2]],tol=1,**kwd,**kwt)
  xa,kwa = run_rf_step(Z0,     tol=.001,x0=[*xd,*xt,*xu])
  log(3,'[x]'.rjust(9)+' ')
  return P

def Rxs_update(P,Pu,q=None,**kwds):
  # update R*x_scen using Pu, e.g. {'Rdx:fsw':1.5}
  # q might be given as quantiles, e.g. Pu[k].ppf(qk) -> 1.5
  if q is not None: Pu = {k:Pu[k].ppf(qk) for k,qk in zip(Pu.keys(),q)}
  for key,sRx in Pu.items():
    pkey,skey = key.split(':')
    s = flatten(strats[skey].ind['s'])
    i = flatten(strats[skey].ind['i'])
    P[pkey+'_scen'][np.ix_(s,i)] = sRx
  P.update(**kwds)
  return P

def run_ss(Ns=10,seed=0):
  log(0,'art.run_ss')
  P0s = fio.load_npy(fname('npy','fit','Ps'))
  Ps = get_sens_sample(P0s,Ns,seed=seed)
  Rs = system.run_n(Ps,t=tvec['main'])
  fio.save_csv(fname('csv','art-ss','wiw',case='sens'),out.wiw(Rs,**tkp))
  fio.save_csv(fname('csv','art-ss','P0s',case='sens'),get_par_expo(P0s,
    keys=[*params.def_sample_distrs().keys(),'PX_fsw','PX_cli','EHY_acute']))
  # TODO: base expo too?
  fio.save_csv(fname('csv','art-ss','expo',case='sens'),merge_expo([
      out.expo([R for R in Rs if R['P']['ss']==ss],ecols=dict(ss=ss),**tko,**ekwds)
      for ss in range(Ns)]))
  fit.plot_sets(tvec['main'],Rs,tfname=fname('fig','art-ss','{}',case='sens'),
    sets='cascade',skeys=['all','aq','fsw','cli'])

def get_par_expo(Ps,keys):
  return dict(par=keys,**{'i'+str(P['id']):[P[key] for key in keys] for P in Ps})

def merge_expo(Es):
  return {k:[x for E in Es for x in E[k]] for k in Es[0]}

def get_sens_sample(Ps,Ns,seed):
  D0 = {
    'Rdx': stats.betabin(p=.65,n=5.3), # CI: (.25,.95)
    'Rtx': stats.betabin(p=.65,n=5.3), # CI: (.25,.95)
    'Rux': stats.gamma(m=6.5,sd=3.5),  # CI: (1.5, 15)
  }
  D = {pkey+':'+skey: D0[pkey] for skey in ('fsw','cli','aq') for pkey in D0}
  return [Rxs_update(deepcopy(P),Pu,ss=ss)
    for ps,P  in enumerate(Ps)
    for ss,Pu in enumerate(params.get_n_sample_lhs(D,Ns,seed=seed+ps))]
