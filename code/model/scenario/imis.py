import numpy as np
from utils import log,fio,stats,parallel,dict_list_update
from model import system,params,target,fit
from model.scenario import N,tvec,fname,get_seeds

PD = params.def_sample_distrs()
PDk = PD.keys()
pic = np.diag([1/D.var() for D in PD.values()]) # prior inverse cov (ignores constr)

# TODO: stopping criterion

def get_ll_wts(lls):
  # hack transform of badly-scaled lls
  return np.quantile(lls,N['fisam']) / np.array(lls) * -np.sign(lls)

def rescale(wts):
  return np.array(wts) / np.sum(wts)

def param_array(Ps):
  return np.array([[P[k] for k in PD.keys()] for P in Ps])

def param_dict(Pa,**kwds):
  return dict_list_update([dict(zip(PD.keys(),Pai)) for Pai in Pa],**kwds)

def get_weights(Rs,Pa,Gs):
  log(2,'imis.get_weights: N = '+str(len(Rs)))
  wp  = N['bsam'] / len(Rs)
  lls = get_ll_wts([R['ll'] for R in Rs]) # likelihood
  lps = np.array([params.get_lp(R['P'],PD) for R in Rs]) # sampling prior
  lgs = np.sum([G.logpdf(Pa) for G in Gs],axis=0) # mvn prior
  lqs = np.log((wp)*np.exp(lps) + (1-wp)*np.exp(lgs)) # mixture prior
  return rescale(lls + lps - lqs) # overall weight

def get_mvn(Pa,wts):
  # get mvn around current best params
  z = np.argmax(wts)
  log(2,'imis.get_gauss: z = '+str(z))
  Pza = Pa[z,] # best params
  dPa = ((Pa-Pza) @ pic * (Pa-Pza)).sum(axis=1) # mahalanobis^2 wrt prior
  zs = np.argsort(dPa)[:N['bisam']+1] # closest params
  Pzcov = np.cov(Pa[zs,].T,aweights=wts[zs]+1/len(wts))
  # assert np.linalg.matrix_rank(Pzcov) == len(PD) # DEBUG
  return stats.mvn(Pza,Pzcov)

def sample_gauss(G,i,gsam=10,jmax=100):
  # sample from G but ensure valid samples (prior > 0)
  log(2,'imis.sample_gauss: N = '+str(N['bisam']))
  # expensive so we do in parallel + batches (gsam)
  def sample_fun(z):
    for j in range(jmax): # attempt
      seed = N['bsam'] + N['bisam']*jmax*i + jmax*z + j
      Ps  = param_dict(G.rvs(gsam,random_state=seed),seed=seed)
      lps = [params.get_lp(P) for P in Ps]
      if np.max(lps) > -np.inf: break # success: prior > 0
    log(3,str(seed).rjust(9)+' ')
    return params.get_all(Ps[np.argmax(lps)])
  return log(-1,parallel.ppool(N['bisam']).map(sample_fun,range(N['bisam'])))

def run(case,b,**kwds):
  log(0,'imis.run: b = {}'.format(b))
  # initialize & first run
  seeds = get_seeds(b)
  T  = target.get_all_esw()
  Gs = []
  Ps = params.get_n_all(seeds,**kwds)
  Rs = system.run_n(Ps,t=tvec['cal'],T=T)
  wts = rescale(get_ll_wts([R['ll'] for R in Rs]))
  # iterations
  for i in range(N['imis']):
    log(1,'imis.iter: i = {}'.format(i))
    zi = slice(len(Ps),len(Ps)+N['bisam'])
    Gs += [get_gauss(param_array(Ps),wts)]
    Ps += sample_gauss(Gs[-1],i)
    Rs += system.run_n(Ps[zi],t=tvec['cal'],T=T)
    wts = get_weights(Rs,param_array(Ps),Gs)
  # saving etc. - TODO: clean-up
  lkwds = dict(wt=wts,ll=[R['ll'] for R in Rs])
  fio.save_csv(fname('csv','imis','Ps',b=b),param_dict(param_array(Ps),lkwds=lkwds))
  Rp = np.random.choice(Rs,size=30,replace=True,p=wts)
  fit.plot_sets(tvec['cal'],Rp,T=T,debug=True)
