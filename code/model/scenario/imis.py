import numpy as np
from utils import log,fio,stats,parallel,dict_list_update
from model import system,params,target,fit
from model.scenario import akwds,N,tvec,fname,get_seeds

# TODO: stopping criterion

def get_ll_wts(lls):
  # hack transform of badly-scaled lls
  # TODO: quantile changes across batches
  return np.quantile(lls,N['isam']/N['bsam']) / np.array(lls) * -np.sign(lls)

def rescale(wts):
  return np.array(wts) / np.sum(wts)

def P_array(Ps,PX):
  return np.array([[P[k] for k in PX['ks']] for P in Ps])

def P_dict(Pa,PX,**kwds):
  return dict_list_update([dict(zip(PX['ks'],Pai)) for Pai in Pa],**kwds)

def init_weights(PX,Rs):
  log(2,'imis.init_weights: N = '+str(len(Rs)))
  W = dict()
  W['ll'] = np.array([R['ll'] for R in Rs])
  W['lp'] = np.array([params.get_lp(R['P'],PX['PD']) for R in Rs])
  W['wt'] = rescale(get_ll_wts(W['ll']))
  return W

def update_weights(W,PX,Rs,zi):
  log(2,'imis.update_weights: N = '+str(len(Rs)))
  wp = N['bsam'] / len(Rs)
  Pa = P_array([R['P'] for R in Rs],PX)
  W['ll'] = np.append(W['ll'],[R['ll'] for R in Rs[zi]]) # likelihood
  W['lp'] = np.append(W['lp'],[params.get_lp(R['P'],PX['PD']) for R in Rs[zi]]) # original prior
  W['lg'] = np.sum([G.logpdf(Pa) for G in PX['GD']],axis=0) # mvn prior
  W['lq'] = np.log((wp)*np.exp(W['lp']) + (1-wp)*np.exp(W['lg'])) # mixture prior
  W['wt'] = rescale(W['ll'] + W['lp'] - W['lq']) # overall weight

def get_mvn(PX,W,Ps):
  # get mvn around current best params
  z = np.argmax(W['wt'])
  log(2,'imis.get_mvn: z = '+str(z))
  Pa = P_array(Ps,PX)
  Pza = Pa[z,] # best params
  dPa = ((Pa-Pza) @ PX['pic'] * (Pa-Pza)).sum(axis=1) # mahalanobis^2 wrt prior
  zs = np.argsort(dPa)[:N['isam']+1] # closest params
  Pzcov = np.cov(Pa[zs,].T,aweights=W['wt'][zs]+1/len(Ps))
  # assert np.linalg.matrix_rank(Pzcov) == len(PX['PD']) # DEBUG
  return stats.mvn(Pza,Pzcov)

def sample_mvn(PX,i,gsam=10,jmax=100):
  # sample from G but ensure valid samples (prior > 0)
  log(2,'imis.sample_mvn: N = '+str(N['isam']))
  # expensive so we do in parallel + batches (gsam)
  def sample_fun(z):
    for j in range(jmax): # attempt
      seed = N['bsam'] + N['isam']*jmax*i + jmax*z + j
      Ps  = P_dict(PX['GD'][i].rvs(gsam,random_state=seed),PX,seed=seed)
      lps = [params.get_lp(P) for P in Ps]
      if np.max(lps) > -np.inf: break # success: prior > 0
    log(3,str(seed).rjust(9)+' ')
    return params.get_all(Ps[np.argmax(lps)])
  return log(-1,parallel.ppool(N['isam']).map(sample_fun,range(N['isam'])))

def init_PX(PD):
  PX = dict()
  PX['PD'] = PD # distributions
  PX['ks'] = list(PD.keys()) # names
  PX['pic'] = np.diag([1/D.var() for D in PD.values()]) # prior inverse cov (ignores constr)
  PX['GD'] = [] # mvn distributions
  return PX

def run(case,b,**kwds):
  log(0,'imis.run: b = {}'.format(b))
  # initialize & first run
  seeds = get_seeds(b)
  PX = init_PX(params.def_sample_distrs())
  T  = target.get_all_esw()
  Ps = params.get_n_all(seeds,**kwds)
  Rs = system.run_n(Ps,t=tvec['cal'],T=T)
  W  = init_weights(PX,Rs)
  # iterations
  for i in range(N['imis']):
    log(1,'imis.iter: i = {}'.format(i))
    zi = slice(len(Ps),len(Ps)+N['isam'])
    PX['GD'] += [get_mvn(PX,W,Ps)]
    Ps += sample_mvn(PX,i)
    Rs += system.run_n(Ps[zi],t=tvec['cal'],T=T)
    update_weights(W,PX,Rs,zi)
  # saving etc. - TODO: clean-up
  lkwds = dict(seed=[P['seed'] for P in Ps],ll=W['ll'],w=W['wt'])
  fio.save_csv(fname('csv','imis','Ps',b=b),P_dict(P_array(Ps,PX),PX,lkwds=lkwds))
  Rp = np.random.choice(Rs,size=100,replace=True,p=W['wt'])
  fit.plot_sets(tvec['cal'],Rp,T=T,debug=True)

if __name__ == '__main__':
  run(**akwds)
  pass
