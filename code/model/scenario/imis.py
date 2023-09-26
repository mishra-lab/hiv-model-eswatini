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

def param_array(Ps,ks):
  return np.array([[P[k] for k in ks] for P in Ps])

def param_dict(Pa,ks,**kwds):
  return dict_list_update([dict(zip(ks,Pai)) for Pai in Pa],**kwds)

def get_weights(Rs,Pa,Gs,PD):
  log(2,'imis.get_weights: N = '+str(len(Rs)))
  wp  = N['bsam'] / len(Rs)
  lls = get_ll_wts([R['ll'] for R in Rs]) # likelihood
  lps = np.array([params.get_lp(R['P'],PD) for R in Rs]) # sampling prior
  lgs = np.sum([G.logpdf(Pa) for G in Gs],axis=0) # mvn prior
  lqs = np.log((wp)*np.exp(lps) + (1-wp)*np.exp(lgs)) # mixture prior
  return rescale(lls + lps - lqs) # overall weight

def get_mvn(Pa,pic,wts):
  # get mvn around current best params
  z = np.argmax(wts)
  log(2,'imis.get_mvn: z = '+str(z))
  Pza = Pa[z,] # best params
  dPa = ((Pa-Pza) @ pic * (Pa-Pza)).sum(axis=1) # mahalanobis^2 wrt prior
  zs = np.argsort(dPa)[:N['isam']+1] # closest params
  Pzcov = np.cov(Pa[zs,].T,aweights=wts[zs]+1/len(wts))
  # assert np.linalg.matrix_rank(Pzcov) == len(PD) # DEBUG
  return stats.mvn(Pza,Pzcov)

def sample_mvn(G,i,ks,gsam=10,jmax=100):
  # sample from G but ensure valid samples (prior > 0)
  log(2,'imis.sample_mvn: N = '+str(N['isam']))
  # expensive so we do in parallel + batches (gsam)
  def sample_fun(z):
    for j in range(jmax): # attempt
      seed = N['bsam'] + N['isam']*jmax*i + jmax*z + j
      Ps  = param_dict(G.rvs(gsam,random_state=seed),ks,seed=seed)
      lps = [params.get_lp(P) for P in Ps]
      if np.max(lps) > -np.inf: break # success: prior > 0
    log(3,str(seed).rjust(9)+' ')
    return params.get_all(Ps[np.argmax(lps)])
  return log(-1,parallel.ppool(N['isam']).map(sample_fun,range(N['isam'])))

def run(case,b,**kwds):
  log(0,'imis.run: b = {}'.format(b))
  # initialize & first run
  seeds = get_seeds(b)
  PD = params.def_sample_distrs()
  ks = list(PD.keys())
  pic = np.diag([1/D.var() for D in PD.values()]) # prior inverse cov (ignores constr)
  T  = target.get_all_esw()
  Gs = []
  Ps = params.get_n_all(seeds,**kwds)
  Rs = system.run_n(Ps,t=tvec['cal'],T=T)
  wts = rescale(get_ll_wts([R['ll'] for R in Rs]))
  # iterations
  for i in range(N['imis']):
    log(1,'imis.iter: i = {}'.format(i))
    zi = slice(len(Ps),len(Ps)+N['isam'])
    Gs += [get_mvn(param_array(Ps,ks),pic,wts)]
    Ps += sample_mvn(Gs[-1],i,ks)
    Rs += system.run_n(Ps[zi],t=tvec['cal'],T=T)
    wts = get_weights(Rs,param_array(Ps,ks),Gs)
  # saving etc. - TODO: clean-up
  lkwds = dict(seed=[P['seed'] for P in Ps],ll=[R['ll'] for R in Rs],wt=wts)
  fio.save_csv(fname('csv','imis','Ps',b=b),param_dict(param_array(Ps,ks),ks,lkwds=lkwds))
  Rp = np.random.choice(Rs,size=100,replace=True,p=wts)
  fit.plot_sets(tvec['cal'],Rp,T=T,debug=True)

if __name__ == '__main__':
  run(**akwds)
  pass
