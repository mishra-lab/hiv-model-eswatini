import re
import numpy as np
from copy import deepcopy
from utils import stats,log,fio,parallel,flatten,minimize
from model import slicers,params,system,target,fit,out
from model.scenario import akwds,N,tvec,fname,batch_select

cascade = dict( # 2020 cascade targets
  low  = (.60,.40,.80),
  mid  = (.80,.80,.90),
  high = (.95,.95,.95),
)
cases = ['fsw+cli+','fsw+cli-','fsw-cli+','fsw-cli-']

def parse_case(case):
  return re.findall('(.*?)(\+|\-)',case)

# TODO: add back keyout

def run_rf(case,b):
  log(0,'scenario.art.run_rf: '+case+' '+str(b))
  P0s = batch_select(fio.load(fname('npy','fit','Ps',case='base')),b=b)
  T = get_refit_T(case+'all-')
  PD = get_refit_PD(case+'aq-')
  fun = lambda P: run_rf_1(P,PD,T,tvec['cal'],ftol=.1)
  Ps = parallel.ppool(len(P0s)).map(fun,P0s); log(1)
  fio.save(fname('npy','art-rf','Ps',case=case,b=b),Ps)

def merge_rf():
  log(0,'scenario.art.merge_rf')
  for case in cases:
    Ps = [P for b in range(N['batch']) for P in fio.load(fname('npy','art-rf','Ps',case=case,b=b))]
    fio.save(fname('npy','art-rf','Ps',case=case),Ps)

def rerun_rf():
  log(0,'scenario.art.rerun_rf')
  tkwds = dict(tvec=tvec['main'],t=tvec['plot'])
  for case in ['base']+cases:
    log(1,case); base = case=='base'
    T = get_refit_T('fsw+cli+all+' if base else case+'all-')
    Ps = fio.load(fname('npy','fit' if base else 'art-rf','Ps',case=case))
    R1s = system.run_n(Ps,t=tvec['main'],T=T)
    fio.save_csv(fname('csv','art-rf','wiw',case=case),out.wiw(R1s,**tkwds))
    fio.save_csv(fname('csv','art-rf','expo',case=case),out.expo(R1s=R1s,**tkwds,
      snames=['all','w','m','aq','fsw','cli'],
      onames=['incidence','prevalence','diagnosed','treated_c','treated_u','vls_c','vls_u']))
    # TODO: add back keyout
    fit.plot_sets(tvec['main'],R1s,T=T,fname=fname('fig','art-rf','{}',case=case),
      sets='cascade',snames=['all','aq','fsw','cli'])
    if base:
      R2s = deepcopy(R1s)
    else:
      fio.save_csv(fname('csv','art-rf','wiw-diff',case=case),out.wiw(R1s,**tkwds,R2s=R2s,vsop='1-2'))

def get_refit_T(case):
  Ts = {
    'fsw-': target.make_targets_2020(cascade['low'], s=0,i=(2,3)),
    'fsw+': target.make_targets_2020(cascade['high'],s=0,i=(2,3),w=1e-6),
    'cli-': target.make_targets_2020(cascade['low'], s=1,i=(2,3)),
    'cli+': target.make_targets_2020(cascade['high'],s=1,i=(2,3),w=1e-6),
    'all-': target.make_targets_2020(cascade['mid'], s=(0,1),i=(0,1,2,3)),
    'all+': target.make_targets_2020(cascade['high'],s=(0,1),i=(0,1,2,3),w=1e-6),
  }
  return flatten( Ts[pop+c] for pop,c in parse_case(case) )

def get_refit_PD(case):
  eps = 1e-9
  PDs = {
    'd+': stats.uniform(l=1-eps,h=1+eps), 'd-': stats.uniform(l=0,h=1),
    't+': stats.uniform(l=1-eps,h=1+eps), 't-': stats.uniform(l=0,h=1),
    'u+': stats.uniform(l=1-eps,h=1+eps), 'u-': stats.uniform(l=1,h=20),
  }
  return {'R'+step+'x:'+pop: PDs[step+c] for pop,c in parse_case(case) for step in 'dtu'}

def run_rf_1(P,PD,T,t,ftol=.1):
  kwds = dict(method='SLSQP',options=dict(ftol=ftol))
  PDz,x0z,Tz,P['refit_iter'] = {},[],[],[]
  for pz,oz in zip(['Rdx','Rtx','Rux'],['diagnosed','treated','vls']):
    PDzi = {k:v for k,v in PD.items() if (pz in k)} # new fitted params
    PDz.update(PDzi)                                # add to fitted params
    x0z += [.5]*len(PDzi)                           # add to initial values
    Tz  += [Ti for Ti in T if (oz in Ti.name)]      # add to targets
    jfun = lambda x: - system.run(Rxs_update(P,PDz,x),t,T,RPts=[])['ll']
    M = minimize(jfun,x0z,bounds=[(0,1) for k in PDz],**kwds)
    x0z[:] = M.x
    if not M.success: return False
  return P

def Rxs_update(P,Pu,q=None):
  # update R*x_scen using Pu, e.g. {'Rdx:fsw':1.5}
  # q might be given as quantiles, e.g. Pu[k].ppf(qk) -> 1.5
  def Rx_update_S(Rx,sname,sRx):
    s = flatten(slicers[sname].pop['s'])
    i = flatten(slicers[sname].pop['i'])
    Rx[np.ix_(s,i)] = sRx
    return Rx
  if q is not None: Pu = {k:Pu[k].ppf(qk) for k,qk in zip(Pu.keys(),q)}
  for key,sRx in Pu.items():
    rate,pop = key.split(':')
    P.update({rate+'_scen':Rx_update_S(P[rate+'_scen'],pop,sRx)})
  return P

def run_ss(Ns=10):
  log(0,'scenario.art.run_ss')
  tkwds = dict(tvec=tvec['main'],t=tvec['plot'])
  P0s = fio.load(fname('npy','fit','Ps'))
  Ps = get_sens_sample(P0s,Ns)
  Rs = system.run_n(Ps,t=tvec['main'])
  fio.save_csv(fname('csv','art-ss','wiw',case='sens'),out.wiw(Rs,**tkwds))
  # TODO: add back keyouts
  fit.plot_sets(tvec['main'],Rs,fname=fname('fig','art-ss','{}',case='sens'),
    sets='cascade',snames=['all','aq','fsw','cli'])

def get_sens_sample(Ps,Ns):
  PDs = { # TODO: adjust to get desire posterior cascade
    'd': stats.betabin(p=.6,n=3),
    't': stats.betabin(p=.6,n=3),
    'u': stats.gamma(m=6,sd=3),
  }
  PD = {'R'+step+'x:'+pop: PDs[step] for pop in ('fsw','cli','aq') for step in 'dtu'}
  return [Rxs_update(deepcopy(P),Pu) for P in Ps
    for Pu in params.get_n_sample_lhs(Ns,PD,seed=P['seed'])]

if __name__ == '__main__':
  # run_rf(**akwds)
  # merge_rf()
  # rerun_rf()
  # run_ss(**akwds)
  pass
