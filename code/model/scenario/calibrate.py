import numpy as np
from utils import log
from model import target,params,system,fit,out
from model.scenario import *

case = 'base'

def run():
  log(0,'scenario.calibrate.run')
  # get targets (T), seeds, params (P0s)
  T = target.get_all_esw()
  seeds = batch_select(range(N['cal']))
  P0s = params.get_n_all(len(seeds),seeds=seeds)
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
  P0s = params.get_n_all(N['cal'],seeds=seeds,all=False)
  for b in range(N['batch']):
    for ll in fio.load_csv(fname('csv','cal','ll',case=case,b=b),fmt='dict'):
      P0s[seeds.index(int(ll['seed']))].update(ll=ll['ll'])
  fio.save_csv(fname('csv','cal','Ps',case=case,b='all'),P0s)

def rerun():
  log(0,'scenario.calibrate.rerun')
  T = target.get_all_esw()
  Ps = fio.load(fname('npy','fit','Ps',case=case,b='all'))
  Rs = system.run_n(Ps,t=tvec['main'],T=T)
  # fit.plot_cal(tvec['main'],Rs,T,fname=fname('fig','fit','plot',case=case,b='all'))
  fio.save_csv(fname('csv','fit','infs',case=case,b='all'),
    out.get_infections(Rs,tvec['main'],tvec['infs']))
