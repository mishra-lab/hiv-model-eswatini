import re
import numpy as np
from copy import deepcopy
from utils import stats,log,fio,parallel,flatten,minimize
from model import strats,params,system,target,fit,out
from model.scenario import akwds,tvec,fname

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

def parse_case(case):
  return re.findall('(.*?)(\+|\-)',case)

def run_rf(b):
  case = cases[b]
  log(0,'art.run_rf: {}'.format(case))
  P0s = fio.load(fname('npy','fit','Ps',case='base'))
  T = get_refit_T(case+'all-')
  D = get_refit_D(case+'aq-')
  fun = lambda P: run_rf_1(P,D,T,tvec['cal'],ftol=.1)
  Ps = parallel.ppool(len(P0s)).map(fun,P0s); log(1)
  fio.save(fname('npy','art-rf','Ps',case=case),Ps)

def rerun_rf():
  log(0,'art.rerun_rf')
  for case in ['base']+cases:
    log(1,case)
    base = (case == 'base')
    T = get_refit_T('fsw+cli+all+' if base else case+'all-')
    Ps = fio.load(fname('npy','fit' if base else 'art-rf','Ps',case=case))
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
  Ts = {
    'fsw-': target.make_targets_2020(cascade['low'], s=0,i=(2,3)),
    'fsw+': target.make_targets_2020(cascade['high'],s=0,i=(2,3),w=1e-6),
    'cli-': target.make_targets_2020(cascade['low'], s=1,i=(2,3)),
    'cli+': target.make_targets_2020(cascade['high'],s=1,i=(2,3),w=1e-6),
    'all-': target.make_targets_2020(cascade['mid'], s=(0,1),i=(0,1,2,3)),
    'all+': target.make_targets_2020(cascade['high'],s=(0,1),i=(0,1,2,3),w=1e-6),
  }
  return flatten( Ts[skey+c] for skey,c in parse_case(case) )

def get_refit_D(case):
  eps = 1e-9
  Ds = {
    'd+': stats.uniform(l=1-eps,h=1+eps), 'd-': stats.uniform(l=0,h=1),
    't+': stats.uniform(l=1-eps,h=1+eps), 't-': stats.uniform(l=0,h=1),
    'u+': stats.uniform(l=1-eps,h=1+eps), 'u-': stats.uniform(l=1,h=20),
  }
  return {'R'+step+'x:'+skey: Ds[step+c] for skey,c in parse_case(case) for step in 'dtu'}

def run_rf_1(P,D,T,t,ftol=.1):
  kwds = dict(method='SLSQP',options=dict(ftol=ftol))
  Dz,x0z,Tz = {},[],[]
  for pz,oz in zip(['Rdx','Rtx','Rux'],['diagnosed','treated','vls']):
    Dzi = {k:v for k,v in D.items() if (pz in k)} # new fitted params
    Dz.update(Dzi)                                # add to fitted params
    x0z += [.5]*len(Dzi)                          # add to initial values
    Tz  += [Ti for Ti in T if (oz in Ti.name)]    # add to targets
    jfun = lambda x: - system.run(Rxs_update(P,Dz,x),t,T,RPts=[])['ll']
    M = minimize(jfun,x0z,bounds=[(0,1) for k in Dz],**kwds)
    x0z[:] = M.x
    if not M.success: return False
  return P

def Rxs_update(P,Pu,q=None,**kwds):
  # update R*x_scen using Pu, e.g. {'Rdx:fsw':1.5}
  # q might be given as quantiles, e.g. Pu[k].ppf(qk) -> 1.5
  def Rx_update_S(Rx,skey,sRx):
    s = flatten(strats[skey].ind['s'])
    i = flatten(strats[skey].ind['i'])
    Rx[np.ix_(s,i)] = sRx
    return Rx
  if q is not None: Pu = {k:Pu[k].ppf(qk) for k,qk in zip(Pu.keys(),q)}
  for key,sRx in Pu.items():
    rate,skey = key.split(':')
    P.update({rate+'_scen':Rx_update_S(P[rate+'_scen'],skey,sRx)},**kwds)
  return P

def run_ss(Ns=10,seed=0):
  log(0,'art.run_ss')
  P0s = fio.load(fname('npy','fit','Ps'))
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
  Ds = {
    'd': stats.betabin(p=.65,n=5.3), # CI: (.25,.95)
    't': stats.betabin(p=.65,n=5.3), # CI: (.25,.95)
    'u': stats.gamma(m=6.5,sd=3.5),  # CI: (1.5, 15)
  }
  D = {'R'+step+'x:'+skey: Ds[step] for skey in ('fsw','cli','aq') for step in 'dtu'}
  return [Rxs_update(deepcopy(P),Pu,ss=ss) for P in Ps
    for ss,Pu in enumerate(params.get_n_sample_lhs(D,Ns,seed=seed))]

if __name__ == '__main__':
  # run_rf(**akwds)
  # rerun_rf()
  # run_ss(**akwds)
  pass

