import numpy as np
from utils import log,fio
from model import target,params,system,fit
from model.scenario import N,tvec,fname,batch_select

case = 'base'

def run():
  log(0,'scenario.calibrate.run')
  # get targets (T), seeds, params (P0s)
  T = target.get_all_esw()
  seeds = batch_select(range(N['cal']))
  P0s = params.get_n_all(len(seeds),seeds=seeds)
  keys = ['seed']+list(params.def_sample_distrs().keys())
  fio.save(fname('npy','sam','Ps',case=case),[{k:P[k] for k in keys} for P in P0s])
  # run & save all ll (nan for fails)
  R0s = system.run_n(P0s,t=tvec['cal'],T=T)
  fio.save_csv(fname('csv','cal','ll',case=case),
    {'seed':seeds,'ll':[R['ll'] if R else np.nan for R in R0s]})
  # select, save, plot topcal %
  R0s = system.drop_fails(R0s)[0]
  Rs = target.top_ll(R0s,top=int(len(seeds)*N['topcal']))
  fio.save(fname('npy','cal','Ps',case=case),[R['P'] for R in Rs])
  fit.plot_cal(tvec['cal'],R0s,T,fname=fname('fig','cal','plot',case=case))

def merge():
  log(0,'scenario.calibrate.merge')
  # npy: load topcal, select topfit, save topfit
  P0s = [P for b in range(N['batch']) for P in fio.load(fname('npy','cal','Ps',case=case,b=b))]
  Ps = target.top_ll(P0s,top=int(N['cal']*N['topfit']))
  fio.save(fname('npy','fit','Ps',case=case,b='all'),Ps)
  # csv: get fitted params, add ll, save
  seeds = range(N['cal'])
  P0s = [P for b in range(N['batch']) for P in fio.load(fname('npy','sam','Ps',case=case,b=b))]
  for b in range(N['batch']):
    for ll in fio.load_csv(fname('csv','cal','ll',case=case,b=b),fmt='dict'):
      P0s[seeds.index(int(ll['seed']))].update(ll=ll['ll'])
  fio.save_csv(fname('csv','sam','Ps',case=case,b='all'),P0s)
