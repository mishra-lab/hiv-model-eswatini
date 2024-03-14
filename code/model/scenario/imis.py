import numpy as np
from utils import log,fio,stats,ppool
from model import system,params,target,fit,out
from model.scenario import N,tvec,fname,get_seeds

# TODO: stopping criterion

D = params.def_sample_distrs()

def xform_ll(lls):
  # hack transform of badly-scaled lls
  # WARNING: assume all lls < 0; quantile changes across batches; nan -> -inf
  lls = np.nan_to_num(lls,nan=-np.inf,neginf=-np.inf)
  return np.quantile(lls,N['isam']/N['hsam']) / lls

def rescale(x):
  return np.array(x) / np.sum(x)

def P_array(Ps):
  return np.array([[P[k] for k in D.keys()] for P in Ps])

def P_dict(Pa):
  return [dict(zip(D.keys(),Pai)) for Pai in Pa]

def update_weights(Ps,Rs,Gs,zi):
  log(2,'imis.update_weights: N = '+str(len(Rs)))
  wp = N['hsam'] / len(Rs)
  Pa = P_array(Ps)
  for P,R in zip(Ps[zi],Rs[zi]):
    P.update(ll=R['ll'],lp=params.get_lp(D,P))
  lls = [P['ll'] for P in Ps] # log-likelihoods
  lps = [P['lp'] for P in Ps] # original log-priors
  lgs = np.sum([G.logpdf(Pa) for G in Gs],axis=0) # mvn log-priors
  lqs = np.maximum(-1e6,np.log((wp)*np.exp(lps) + (1-wp)*np.exp(lgs))) # mixture log-priors
  return rescale(xform_ll(lls) * np.exp(lps - lqs)) # overal weight

def get_mvn(wts,Pa):
  # get mvn around current best params
  z = np.argmax(wts)
  log(2,'imis.get_mvn: z = '+str(z))
  Pza = Pa[z,] # best params
  pic = np.diag([1/Dk.var() for Dk in D.values()]) # prior inv cov (ignore constr)
  dPa = ((Pa-Pza) @ pic * (Pa-Pza)).sum(axis=1) # mahalanobis^2 wrt prior
  zs = np.argsort(dPa)[:N['isam']+1] # closest params to best
  Pzcov = np.cov(Pa[zs,].T,aweights=wts[zs]+1/Pa.shape[0])
  # assert np.linalg.matrix_rank(Pzcov) == len(D) # DEBUG
  return [stats.mvn(Pza,Pzcov)]

def sample_mvn(G,gsam=10,jmax=100,**kwds):
  # sample from G but ensure valid samples (prior > 0)
  log(2,'imis.sample_mvn: N = '+str(N['isam']))
  # expensive so we do in parallel + batches (gsam)
  def sample_fun(z):
    for j in range(jmax): # attempt
      Ps  = P_dict(G.rvs(gsam,random_state=jmax*z+j))
      lps = [params.get_lp(D,P) for P in Ps]
      if np.max(lps) > -np.inf: break # success: prior > 0
    log(3,str(z).rjust(9)+' ')
    return params.get_depend(Ps[np.argmax(lps)],id=z,**kwds)
  return log(-1,ppool(N['isam']).map(sample_fun,range(N['isam'])))

def run(case,b,**kwds):
  log(0,'imis.run: {} @ b = {}'.format(case,b))
  # initialize & first run
  seeds = get_seeds(b)
  T = target.get_all_esw()
  Gs = []
  Ps = params.get_n_all(seeds,batch=b,imis=0,**kwds)
  Rs = system.run_n(Ps,t=tvec['cal'],T=T)
  wts = update_weights(Ps,Rs,Gs,slice(N['hsam']))
  # iterations
  for i in range(N['imis']):
    log(1,'imis.iter: i = {}'.format(i+1))
    zi = slice(len(Ps),len(Ps)+N['isam'])
    Gs += get_mvn(wts,P_array(Ps))
    Ps += sample_mvn(Gs[-1],batch=b,imis=i+1,**kwds)
    Rs += system.run_n(Ps[zi],t=tvec['cal'],T=T)
    wts = update_weights(Ps,Rs,Gs,zi)
  kxs = ('id','batch','imis',*D.keys(),'foi_mode','ll','lp')
  Pxs = [dict({k:P[k] for k in kxs},wt=wt) for P,wt in zip(Ps,wts)]
  fio.save_npy(fname('npy','imis','Ps',case=case,b=b),Pxs)
  fio.save_csv(fname('csv','imis','Ps',case=case,b=b),Pxs)

def sample_post(case,seed=0):
  log(0,'imis.sample_post: {}'.format(case))
  P0xs = [P for b in range(N['batch']) for P in fio.load_npy(fname('npy','imis','Ps',case=case,b=b))]
  np.random.seed(seed)
  wll = rescale(xform_ll([P['ll'] for P in P0xs])) # weights
  Pxs = np.random.choice(P0xs,N['post'],replace=False,p=wll) # sample
  for P in Pxs: P.update(id='{}.{}.{}'.format(P['batch'],P['imis'],P['id']))
  fio.save_csv(fname('csv','fit','Ps',case=case),Pxs)
  fio.save_npy(fname('npy','fit','Ps',case=case),[params.get_depend(P) for P in Pxs])

def rerun(case):
  log(0,'imis.rerun: {}'.format(case))
  T = target.get_all_esw()
  Ps = fio.load_npy(fname('npy','fit','Ps',case=case))
  Rs = system.run_n(Ps,t=tvec['main'],Xk=True)
  fio.save_csv(fname('csv','fit','wiw',case=case),out.wiw(Rs,tvec['main'],tvec['plot']))
  fit.plot_sets(tvec['main'],Rs,T=T,tfname=fname('fig','fit','{}',case=case))
