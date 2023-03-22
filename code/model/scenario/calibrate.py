import numpy as np
from utils import log,fio
from model import system,params,target,fit,out
from model.scenario import akwds,N,tvec,fname,batch_select

def run(case,b,**kwds):
  log(0,'scenario.calibrate.run: '+str(b))
  # get targets (T), seeds, params (P0s)
  T = target.get_all_esw()
  seeds = batch_select(range(N['sam']),b=b)
  P0s = params.get_n_all(len(seeds),seeds=seeds,**kwds)
  keys = ['seed','foi_mode']+list(params.def_sample_distrs().keys())
  fio.save(fname('npy','sam','Ps',case=case,b=b),[{k:P[k] for k in keys} for P in P0s])
  # run & save all ll (nan for fails)
  R0s = system.run_n(P0s,t=tvec['cal'],T=T)
  fio.save_csv(fname('csv','sam','ll',case=case,b=b),
    {'seed':seeds,'ll':[R['ll'] if R else np.nan for R in R0s]})
  # select, save, plot topcal %
  R0s = system.drop_fails(R0s)[0]
  Rs = target.top_ll(R0s,top=int(len(seeds)*N['topcal']))
  fio.save(fname('npy','cal','Ps',case=case,b=b),[R['P'] for R in Rs])
  fit.plot_sets(tvec['cal'],R0s,T=T,tfname=fname('fig','sam','cal',case=case,b=b),debug=True)

def merge(case):
  log(0,'scenario.calibrate.merge')
  # npy: load topcal, select topfit, save topfit
  P0s = [P for b in range(N['batch']) for P in fio.load(fname('npy','cal','Ps',case=case,b=b))]
  Ps = target.top_ll(P0s,top=int(N['sam']*N['topfit']))
  fio.save(fname('npy','fit','Ps',case=case),Ps)
  # csv: get fitted params, add ll, save
  seeds = range(N['sam'])
  P0s = [P for b in range(N['batch']) for P in fio.load(fname('npy','sam','Ps',case=case,b=b))]
  for b in range(N['batch']):
    for ll in fio.load_csv(fname('csv','sam','ll',case=case,b=b),fmt='dict'):
      P0s[seeds.index(int(ll['seed']))].update(ll=ll['ll'])
  fio.save_csv(fname('csv','sam','Ps',case=case),P0s)

def rerun(case):
  log(0,'scenario.calibrate.rerun')
  T = target.get_all_esw()
  Ps = fio.load(fname('npy','fit','Ps',case=case))
  Rs = system.run_n(Ps,t=tvec['main'])
  fit.plot_sets(tvec['main'],Rs,T=T,tfname=fname('fig','fit','{}',case=case))
  fio.save_csv(fname('csv','fit','wiw',case=case),out.wiw(Rs,tvec['main'],tvec['plot']))

if __name__ == '__main__':
  # run(**akwds)
  # merge(**akwds)
  # rerun(**akwds)
  pass
