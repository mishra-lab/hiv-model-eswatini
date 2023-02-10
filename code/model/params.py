import numpy as np
from copy import copy
from scipy.optimize import nnls
from utils import _,NAN,log,stats,flatten,dict_split,linear_comb,interval_qs
from utils import tarray as ta
from model import tol

# main -------------------------------------------------------------------------

def get_all(P=None,seed=None,**kwds):
  P = P if P is not None else \
      get_sample_random(seed=seed)
  P['foi_mode'] = 'base'
  P.update(kwds)
  # order matters for some dependencies
  P.update(get_F(P))
  P.update(get_PX(P))
  P.update(get_birth_death(P))
  # P.update(get_turnover_old(P))
  P.update(get_turnover(P))
  P.update(get_K(P))
  P.update(get_mix(P))
  P.update(get_beta_a(P))
  P.update(get_condom(P))
  P.update(get_circumcision(P))
  P.update(get_hiv_prog(P))
  P.update(get_diag(P))
  P.update(get_treat(P))
  P.update(get_scen(P))
  return P

def get_n_all(n,Ps=None,seeds=None,lhs=True,all=True,**kwds):
  log(2,'params.get_n_all: '+str(n))
  if seeds is None: seeds = n*[None]
  if Ps is None:
    PD = def_sample_distrs()
    if lhs:
      # only latin hypercube sample params without constraints
      # too expensive / frail for constraints - https://arxiv.org/abs/0909.0329
      PDc = dict_split(PD,flatten(def_checkers().values()))
      Phs = get_n_sample_lhs(n,PD,seed=seeds[0])
      Ps = get_n_sample_random(n,PDc,seeds=seeds,Ps=Phs)
    else:
      Ps = get_n_sample_random(n,PD,seeds=seeds)
  if all:
    return [get_all(P,**kwds) for P in Ps]
  else:
    return Ps

def get_n_sample_lhs(n,PD,seed=None):
  Qs = stats.lhs(len(PD),n,seed=seed)
  return [{key:PD[key].ppf(q) for key,q in zip(PD,Q)} for Q in Qs]

def get_n_sample_random(n,PD=None,seeds=None,Ps=None):
  if Ps is None: Ps = n*[{}]
  if seeds is None: seeds = n*[None]
  if PD is None: PD = def_sample_distrs()
  return [get_sample_random(PD,seed=seeds[i],P=Ps[i]) for i in range(n)]

def get_sample_random(PD=None,seed=None,P=None):
  if P is None: P = {}
  if seed is not None: np.random.seed(seed)
  if PD is None: PD = def_sample_distrs()
  P.update({key:dist.rvs() for key,dist in PD.items()})
  P['seed'] = seed
  checkers = def_checkers()
  # adjustments / forcing
  for checker,keys in checkers.items():
    resample_until(P,PD,checker,keys)
  return P

def resample_until(P,PD,checker,keys,log=False):
  if keys is None: keys = PD.keys()
  n = 1
  while not checker(P):
    n += 1
    for key in keys:
      P[key] = PD[key].rvs()
  if log:
    print('{:.5f} {:>15} @ {}'.format((1/n),str(checker).split()[1],', '.join(keys)))
  return P

def def_checkers():
  k = 'PF_condom_' # convenience
  return {
    check_F:      ['PF_ai_mcx','PF_ai_swx','F_swr'],
    check_acute:  ['Rbeta_acute','dur_acute'],
    check_gud:    ['P_gud_fsw_l','RP_gud_fsw_h:l'],
    check_condom: [k+'msp_2006',k+'msp_2016',
                   k+'cas_2006',k+'cas_2016',
                   k+'swo_2002',k+'swo_2011',k+'swo_2014',
                   k+'swr_2002',k+'swr_2011',k+'swr_2014']
  }

def def_sample_distrs():
  return {
  # pop
  't0_hiv':               stats.uniform(l=1980,h=1985),
  'PX_w_fsw':             stats.betabin(p=.0288,n=121),
  'PX_w_h':               stats.betabin(p=.1778,n=66),
  'PX_m_m':               stats.betabin(p=.1331,n=360),
  'dur_fsw_l':            stats.gamma(m=4.125,sd=0.644),
  'dur_fsw_h':            stats.gamma(m=9.408,sd=1.748),
  'dur_cli':              stats.gamma(m=8.626,sd=2.832),
  'turn_xm_xl':           stats.gamma(m=.2156,sd=.1177),
  'Pturn_fsw_m:l':        stats.betabin(p=.724,n=18),
  'Pturn_cli_m:l':        stats.betabin(p=.602,n=7),
  'growth_2050':          stats.uniform(l=.007,h=.015),
  # partners
  'C12m_msp_xl':          stats.betabin(p=.370,n=57),
  'C12m_cas_xl':          stats.betabin(p=.366,n=28),
  'C12m_cas_wm':          stats.gamma(m=1.575,sd=0.204),
  'RC_cas_cli:wm':        stats.uniform(l=.25,h=1),
  'C1m_swo_fsw_l':        stats.gamma(m=4.06,sd=0.897),
  'C1m_swr_fsw_l':        stats.gamma(m=6.93,sd=1.065),
  'RC_swo_fsw_h:l':       stats.gamma(m=2.02,sd=0.230),
  'RC_swr_fsw_h:l':       stats.gamma(m=1.49,sd=0.102),
  'KF_swx_cli':           stats.gamma(m=59.3,sd=14.1),
  'RKF_swx_cli_h:l':      stats.gamma(m=2.03,sd=0.230),
  # sex
  'F_msp':                stats.gamma(m=77.3,sd=33.7),
  'RF_cas:msp':           stats.uniform(l=.25,h=1),
  'dur_msp':              stats.uniform(l=14.5,h=18.5),
  'dur_cas':              stats.gamma(m=.743,sd=0.324),
  'dur_swr':              stats.gamma(m=1.125,sd=0.386),
  'F_swr':                stats.uniform(l=12,h=36),
  'PF_ai_mcx':            stats.gamma(m=.0573,sd=.0423),
  'PF_ai_swx':            stats.gamma(m=.0920,sd=.0769),
  # mixing
  'pref_msp_xl':          stats.gamma(m=2.19,sd=.384),
  'pref_mcx_swx':         stats.gamma(m=8.3,sd=4.44),
  # condoms
  'Rbeta_condom':         stats.betabin(p=.734,n= 32),
  'RPF_condom_a:v':       stats.betabin(p=.768,n= 12),
  'PF_condom_msp_2006':   stats.betabin(p=.230,n=100),
  'PF_condom_msp_2016':   stats.betabin(p=.416,n= 75),
  'PF_condom_cas_1988':   stats.betabin(p=.042,n= 40),
  'PF_condom_cas_2006':   stats.betabin(p=.598,n=100),
  'PF_condom_cas_2016':   stats.betabin(p=.694,n=100),
  'PF_condom_swo_2002':   stats.betabin(p=.432,n=  9),
  'PF_condom_swo_2011':   stats.betabin(p=.777,n= 21),
  'PF_condom_swo_2014':   stats.betabin(p=.787,n= 14),
  'PF_condom_swr_2002':   stats.betabin(p=.337,n= 13),
  'PF_condom_swr_2011':   stats.betabin(p=.754,n= 24),
  'PF_condom_swr_2014':   stats.betabin(p=.759,n= 11),
  'PF_circum_2050':       stats.betabin(p=.724,n= 18),
  # beta
  'beta_0':               stats.gamma(m=.00131,sd=.00052),
  'Rbeta_acute':          stats.gamma(m=6.015,sd=3.653),
  'Rbeta_350':            stats.gamma(m=1.586,sd=0.1532),
  'Rbeta_200':            stats.gamma(m=8.203,sd=2.182),
  'Rbeta_vi_rec':         stats.uniform(l=1,h=2),
  'aRbeta_gud_sus':       stats.gamma(m=2.05,sd=1.546),
  'aRbeta_gud_inf':       stats.gamma(m=0.99,sd=0.577),
  'dur_acute':            stats.gamma(m=.172,sd=.1285),
  'P_gud_fsw_l':          stats.betabin(p=.232,n=29),
  'RP_gud_fsw_h:l':       stats.gamma(m=3.00,sd=0.900),
  'RP_gud_2050':          stats.uniform(l=.2,h=1),
  'iP_gud_cli_fsw:gp':    stats.uniform(l=0,h=1),
  'Rbeta_uvls':           stats.betabin(p=.244,n=5),
  # diagnosis
  'Rdx_global':           stats.uniform(l=.5,h=1),
  'dx_w_2002':            stats.betabin(p=.094,n=100),
  'dx_w_2006':            stats.betabin(p=.248,n=100),
  'Rdx_m:w_2006':         stats.gamma(m=.377,sd=.1),
  'dx_wq_2011':           stats.gamma(m=.637,sd=0.094),
  'Rdx_m:wq_2011':        stats.gamma(m=.529,sd=.1),
  'aRdx_fsw:wq_2011':     stats.gamma(m=.5207,sd=.2),
  'aRdx_wq_16:11':        stats.gamma(m=.2035,sd=.05),
  'aRdx_fsw:wq_2016':     stats.gamma(m=.619,sd=.2),
  # treatment
  'tx_2010':              stats.gamma(m=1.5,sd=0.65),
  'tx_2012':              stats.gamma(m=8.75,sd=1.53),
  'Rtx_fsw:wq':           stats.uniform(l=.5,h=1),
  'ivx':                  stats.gamma(m=.62,sd=.172),
  'Runvx_m:wq':           stats.uniform(l=1,h=1.5),
  'Runvx_fsw:wq':         stats.uniform(l=1,h=1.5),
  'revx_2010':            stats.gamma(m=.7288,sd=.1279),
  }

def print_sample_distrs(PD=None,fmt='{:6.3f}',interval=.95,keys=None):
  if PD is None: PD = def_sample_distrs()
  if keys is None: keys = PD.keys()
  r = max(map(len,PD.keys()))
  sf = '{}: '+fmt+' & ('+fmt+',~'+fmt+') | {}'
  for key in keys:
    print(sf.format(key.rjust(r),PD[key].mean(),*PD[key].interval(interval),PD[key].dist.name))

def print_sampled_distrs(Ps,fmt='{:6.3f}',interval=.95,keys=None):
  # TODO: support slicing?
  if keys is None: keys = Ps[0].keys()
  r = max(map(len,keys))
  sf = '{}: '+fmt+' & ('+fmt+',~'+fmt+')'
  for key in keys:
    Pk = [P[key] for P in Ps]
    print(sf.format(key.rjust(r),np.mean(Pk),*np.quantile(Pk,interval_qs(interval))))

# demographics -------------------------------------------------------------------------------------

def get_PX(P): # [OK]
  # population size
  # PX.shape = (s:2, i:4, k:5, h:6, c:5)
  NX0   = np.array(243) # REF: WorldBank
  PX0_k = np.array([1,0,0,0,0]).reshape([1,1,5,1,1])
  PX0_h = np.array([1,0,0,0,0,0]).reshape([1,1,1,6,1])
  PX0_c  = np.array([1,0,0,0,0]).reshape([1,1,1,1,5])
  PX_h_hiv = np.array([0,5,65,30,0,0])*1e-6 # REF: assume
  PX_h_hiv[0] = 1 - PX_h_hiv.sum()
  P = get_PX_swx(P)
  PX_si = np.zeros((2,4))
  # FSW
  PX_si[0,2] = P['PX_fsw'] * (1-P['PX_fsw_h'])
  PX_si[0,3] = P['PX_fsw'] * P['PX_fsw_h']
  # Clients
  PX_si[1,3] = P['PX_cli'] * P['PX_cli_h']
  PX_si[1,2] = P['PX_cli'] * (1-P['PX_cli_h'])
  # non-client men
  PX_si[1,1] = P['PX_m_m'] * (1-P['PX_w'])
  PX_si[1,0] = 1 - P['PX_w'] - PX_si[1,1:].sum()
  # non-FSW women
  PX_si[0,1] = P['PX_w'] * P['PX_w_h'] - P['PX_fsw']
  PX_si[0,0] = P['PX_w'] - PX_si[0,1:].sum()
  return {
    'PX_m_cli': P['PX_cli'] / (1 - P['PX_w']),
    'PX_s':  PX_si.sum(axis=1),
    'PX_si': PX_si,
    'PX_si_s': PX_si / PX_si.sum(axis=1)[:,_],
    'X0': NX0 * PX_si[:,:,_,_,_] * PX0_k * PX0_h * PX0_c,
    'PX_h_hiv': PX_h_hiv,
  }

def get_PX_swx(P):
  P['PX_w']     = .52 # REF: WorldBank
  P['PX_fsw_h'] = .2  # REF: assume
  P['PX_cli_h'] = .2  # REF: assume
  P['F_swo']    = 12  # see get_F()
  P['PX_fsw']   = P['PX_w'] * P['PX_w_fsw']
  P['XKF_swo']  = P['PX_fsw'] * P['C2K_p'][2] * P['C1m_swo_fsw_l'] * P['F_swo'] * \
    linear_comb(P['PX_fsw_h'],P['RC_swo_fsw_h:l'],1)
  P['XKF_swr']  = P['PX_fsw'] * P['C2K_p'][3] * P['C1m_swr_fsw_l'] * P['F_swr'] * \
    linear_comb(P['PX_fsw_h'],P['RC_swr_fsw_h:l'],1)
  P['PX_cli'] = (P['XKF_swo'] + P['XKF_swr']) / P['KF_swx_cli']
  # NOTE: KF_swx_cli is true client average, not KF_swx_cli_l; see get_C()
  return P

def get_birth_death(P): # [OK]
  death = np.array(1/35 + (1-.64)*.0144) # death rate ~ .034 for 15-49 years + non-HIV mortality
  birth = ta.tarray([1980,2000,2010,2020,2050],np.array([.04,.03,.015,.015,P['growth_2050']])+death)
  return {
    'birth_t': birth,
    'death': death,
  }

def get_turnover(P):
  # turn_sii.shape = (s:2, i:4, i':4)
  dur = np.array([[NAN,NAN,P['dur_fsw_l'],P['dur_fsw_h']],[NAN,NAN,P['dur_cli'],P['dur_cli']]])
  PXe = np.zeros((2,4)) # s, i
  turn = np.zeros((2,4,4)) # s, i.from, i.to
  b4 = np.eye(4,dtype=bool) # pre-compute
  A = np.zeros((2,16,16)) #  0, 1, 2, 3,   4,  5,  6,   7,  8,  9,  10, 11, 12,  13, 14, 15
  b = np.zeros((2,16))    # e1,e2,e3,e4, t01,t02,t03, t10,t12,t13, t20,t21,t23, t30,t31,t32
  for s in (0,1):
    x = P['PX_si'][s,:]
    b[s,0:4] = x * NAN # * birth_t (runtime)
    A[s,0,:] = (NAN,0,0,0,-x[0],-x[0],-x[0],+x[1],    0,    0,+x[2],    0,    0,+x[3],    0,    0)
    A[s,1,:] = (0,NAN,0,0,+x[0],    0,    0,-x[1],-x[1],-x[1],    0,+x[2],    0,    0,+x[3],    0)
    A[s,2,:] = (0,0,NAN,0,    0,+x[0],    0,    0,+x[1],    0,-x[2],-x[2],-x[2],    0,    0,+x[3])
    A[s,3,:] = (0,0,0,NAN,    0,    0,+x[0],    0,    0,+x[1],    0,    0,+x[2],-x[3],-x[3],-x[3])
    A[s,4,0:4],b[s,4]   = 1, x.sum() # e.sum = x.sum
    A[s,5,10:13],b[s,5] = 1, 1/dur[s,2] - P['death'] # dur 2
    A[s,6,13:16],b[s,6] = 1, 1/dur[s,3] - P['death'] # dur 3
    A[s,7,12],b[s,7] = 1, 0.0 # t23 = 0
    A[s,8,15],b[s,8] = 1, 0.0 # t32 = 0
    A[s,9, 7],b[s,9] = 1, P['turn_xm_xl'] # t10
    A[s,11,1],b[s,11] = 1, x[1]*1.0 # e1 = x1*R
    A[s,12,2],b[s,12] = 1, NAN # e2 = x2*R (runtime)
    A[s,13,3],b[s,13] = 1, NAN # e3 = x3*R (runtime)
    p = (P['Pturn_fsw_m:l'],P['Pturn_cli_m:l'])[s]
    A[s,14,:] = (0,0,0,0, 0,0,0, 0,0,0, p-1,p,0, 0,0,0) # (p-1) * t20 = (p) * t21
    A[s,15,:] = (0,0,0,0, 0,0,0, 0,0,0, 0,0,0, p-1,p,0) # (p-1) * t30 = (p) * t31
    # print(np.array(sympy.Matrix(A[s,:,:]).rref()[0],dtype=float).round(3)) # DEBUG
  # print(np.linalg.matrix_rank(A)) # DEBUG ~= [16,16]
  return {
    'dur_si': dur,
    'PXe_si': PXe,
    'turn_sii': turn,
    'turn_A_s': A,
    'turn_b_s': b,
    'b4': b4,
  }

def solve_turnover(P,t):
  v = P['birth_t'](t)
  P['turn_b_s'][:,0:4]     = v * P['PX_si']
  P['turn_A_s'][:,0:4,0:4] = v * P['b4']
  for s in (0,1):
    Re = np.minimum((2.0,1.5)[s], (v - P['death'] + 1/P['dur_si'][s,2:4]) / v)
    P['turn_b_s'][s,12:14] = P['PX_si'][s,2:4] * Re
    et,err = nnls(P['turn_A_s'][s],P['turn_b_s'][s])
    P['PXe_si'][s] = et[0:4]
    P['turn_sii'][s,np.logical_not(P['b4'])] = et[4:]
    if err > tol:
      P['PXe_si'][s] = -np.inf # fail (gracefully) in system.solve
      # raise Exception('Cannot solve turnover (s={}) error: {}, seed: {}'.format(s,err,P['seed']))
  return v,P['PXe_si'],P['turn_sii']

def get_turnover_old(P): # TODO: remove
  # turn_sii.shape = (s:2, i:4, i':4)
  t_fsw_l = 1/P['dur_fsw_l'] - P['death']
  t_fsw_h = 1/P['dur_fsw_h'] - P['death']
  t_cli   = 1/P['dur_cli'] - P['death']
  turn = np.zeros((2,4,4)) # s, i.from, i.to
  # fsw -> med & low
  turn[0,2,1] = t_fsw_l * (P['Pturn_fsw_m:l'])
  turn[0,2,0] = t_fsw_l * (1 - P['Pturn_fsw_m:l'])
  turn[0,3,1] = t_fsw_h * (P['Pturn_fsw_m:l'])
  turn[0,3,0] = t_fsw_h * (1 - P['Pturn_fsw_m:l'])
  # clients -> med & low
  turn[1,2,1] = t_cli * (P['Pturn_cli_m:l'])
  turn[1,2,0] = t_cli * (1 - P['Pturn_cli_m:l'])
  turn[1,3,1] = t_cli * (P['Pturn_cli_m:l'])
  turn[1,3,0] = t_cli * (1 - P['Pturn_cli_m:l'])
  # med -> low
  turn[0,1,0] = P['turn_xm_xl']
  turn[1,1,0] = P['turn_xm_xl']
  # balance flows (absolute)
  i = [(3,3,2,2,1),(1,0,1,0,0)]
  turn[:,i[1],i[0]] += turn[:,i[0],i[1]] * P['PX_si'][:,i[0]] / P['PX_si'][:,i[1]]
  return {
    'turn_sii': turn,
    'dur_si': 1 / (turn.sum(axis=2) + P['death']),
  }

# FOI ----------------------------------------------------------------------------------------------

def get_beta_a(P): # [OK]
  # beta_a.shape = (a:2, p:4, s:2, i:4, s':2, i':4, h':6, c':5, (t))
  Rbeta_ar = 10
  Rbeta_as = np.array([[P['Rbeta_vi_rec'],1],[Rbeta_ar,1]]).reshape([2,1,2,1,1,1,1,1])
  Rbeta_as = Rbeta_as / Rbeta_as[0,:].mean()
  P_gud_gp = .07 # REF: SDHS2006 (adj % sex active below)
  P_gud_cli_h = linear_comb(P['iP_gud_cli_fsw:gp'],  P['P_gud_fsw_l'],P_gud_gp)
  P_gud_cli_l = linear_comb(P['iP_gud_cli_fsw:gp']/2,P['P_gud_fsw_l'],P_gud_gp)
  P_gud = np.array([
    [P_gud_gp*.81, P_gud_gp, P['P_gud_fsw_l'], P['P_gud_fsw_l']*P['RP_gud_fsw_h:l']],
    [P_gud_gp*.68, P_gud_gp, P_gud_cli_l, P_gud_cli_h],
  ])
  RP_gud_t = ta.tarray([1980,2000,2020,2050,2051],[1,1,1,*2*[P['RP_gud_2050']]])
  Rbeta_h = np.array([0,P['Rbeta_acute'],1,1,P['Rbeta_350'],P['Rbeta_200']]).reshape([1,1,1,1,1,1,6,1])
  Rbeta_c = np.array([1,1,1-(1-P['Rbeta_uvls'])/2,P['Rbeta_uvls'],.00]).reshape([1,1,1,1,1,1,1,5])
  return {
    'beta_a': P['beta_0'] * Rbeta_as * Rbeta_h * Rbeta_c,
    'Rbeta_as': Rbeta_as,
    'Rbeta_h': Rbeta_h,
    'P_gud': P_gud,
    'RP_gud_t': RP_gud_t,
    'EHY_acute': P['Rbeta_acute'] * P['dur_acute'],
  }

def check_acute(P):
  return 1 <= P['Rbeta_acute'] * P['dur_acute'] <= 63

def check_gud(P):
  return (
    P['P_gud_fsw_l'] > .07 and
    P['P_gud_fsw_l'] * P['RP_gud_fsw_h:l'] < 1
  )

def get_F(P): # [OK]
  # F_ap.shape = (a:2, p:4)
  F_p   = np.array([ P['F_msp'], P['F_msp']*P['RF_cas:msp'], 12, P['F_swr'] ])
  PF_ai = np.array([ P['PF_ai_mcx'],P['PF_ai_mcx'],P['PF_ai_swx'],P['PF_ai_swx'] ])
  F_ap  = F_p.reshape([1,4]) * np.array([1-PF_ai,PF_ai])
  # pcr_msp # TODO: thesis update
  dur_p = np.array([ P['dur_msp'], P['dur_cas'], 1/12, P['dur_swr'] ])
  dur_p_1 = np.minimum(dur_p,1)
  C2K_p = dur_p / (dur_p + np.array([1, 1, 1/12, 1/12]))
  return {
    'F_ap': F_ap,
    'dur_p': dur_p,
    'dur_p_1': dur_p_1,
    'C2K_p': C2K_p,
  }

def check_F(P):
  C2K_swo = 1/2
  C2K_swr = P['dur_swr'] / (P['dur_swr'] + 1/12)
  return (
    P['C1m_swo_fsw_l'] * C2K_swo * P['RC_swo_fsw_h:l'] * 12 + \
    P['C1m_swr_fsw_l'] * C2K_swr * P['RC_swr_fsw_h:l'] * P['F_swr'] < 2*365 and
    P['PF_ai_swx'] <= 0.5 and
    P['PF_ai_mcx'] <= P['PF_ai_swx']
  )

def get_K(P): # TODO C -> Q
  # .shape = (p:4, s:2, i:4)
  PX_si = P['PX_si']
  F = np.squeeze(P['F_ap'].sum(axis=0))
  # dimensions: p,s,i
  K_psi = np.zeros((4,2,4,1))
  # main / spousal
  K_psi[0,0,0] = P['C12m_msp_xl']
  K_psi[0,0,1] = P['C12m_msp_xl']
  K_psi[0,0,2] = 0.5
  K_psi[0,0,3] = 0.5
  K_psi[0,1,0] = P['C12m_msp_xl']
  K_psi[0,1,2] = P['C12m_msp_xl'] * .5
  K_psi[0,1,3] = P['C12m_msp_xl'] * .5
  K_psi[0,1,1] = (K_psi[0,0,:] * PX_si[0,:,_] - K_psi[0,1,:] * PX_si[1,:,_]).sum() / PX_si[1,1]
  K_psi[0] *= P['C2K_p'][0]
  # casual
  K_psi[1,0,0] = P['C12m_cas_xl']
  K_psi[1,0,1] = P['C12m_cas_wm']
  K_psi[1,0,2] = 0.5
  K_psi[1,0,3] = 1.0
  K_psi[1,1,0] = P['C12m_cas_xl']
  K_psi[1,1,2] = P['C12m_cas_wm'] * P['RC_cas_cli:wm']
  K_psi[1,1,3] = P['C12m_cas_wm'] * P['RC_cas_cli:wm']
  K_psi[1,1,1] = (K_psi[1,0,:] * PX_si[1,:,_] - K_psi[1,1,:] * PX_si[1,:,_]).sum() / PX_si[1,1]
  K_psi[1] *= P['C2K_p'][1]
  # swx: fsw
  K_psi[2,0,2] = P['C2K_p'][2] * P['C1m_swo_fsw_l']
  K_psi[2,0,3] = P['C2K_p'][2] * P['C1m_swo_fsw_l'] * P['RC_swo_fsw_h:l']
  K_psi[3,0,2] = P['C2K_p'][3] * P['C1m_swr_fsw_l']
  K_psi[3,0,3] = P['C2K_p'][3] * P['C1m_swr_fsw_l'] * P['RC_swr_fsw_h:l']
  # swx: clients
  wPX = (PX_si[1,2] + P['RKF_swx_cli_h:l'] * PX_si[1,3])
  K_psi[2,1,2] = P['XKF_swo'] / F[2] / wPX
  K_psi[2,1,3] = P['XKF_swo'] / F[2] / wPX * P['RKF_swx_cli_h:l']
  K_psi[3,1,2] = P['XKF_swr'] / F[3] / wPX
  K_psi[3,1,3] = P['XKF_swr'] / F[3] / wPX * P['RKF_swx_cli_h:l']
  KF_psi = K_psi.sum(axis=3) * F[:,_,_]
  aK_pk = np.array([[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]]).reshape((4,1,1,5))
  aK_pk = aK_pk * (K_psi > 0)
  return {
    'K_psi': K_psi,
    'KF_psi': KF_psi,
    'aK_pk': aK_pk,
  }

def get_condom(P): # [OK]
  k = 'PF_condom_' # convenience
  PF_condom_t = ta.tarray([1980,1988,2002,2006,2011,2014,2016,2050],
    [[0,  0,NAN,P[k+'msp_2006'],NAN,NAN,*2*[P[k+'msp_2016']] ], # main
     [0,  0,NAN,P[k+'cas_2006'],NAN,NAN,*2*[P[k+'cas_2016']] ], # casual
     [0,NAN,P[k+'swo_2002'],NAN,P[k+'swo_2011'],*3*[P[k+'swo_2014']] ], # sw-new
     [0,NAN,P[k+'swr_2002'],NAN,P[k+'swr_2011'],*3*[P[k+'swr_2014']] ]] # sw-reg
  ).reshape([1,4,1,1,1,1,1,1])
  RPF_condom_a = np.array([1,P['RPF_condom_a:v']]).reshape([2,1,1,1,1,1,1,1])
  return {
    'PF_condom_t': PF_condom_t,
    'RPF_condom_a': RPF_condom_a,
  }

def check_condom(P):
  # TODO: NERCHA2018: cas condom use may be declining?
  k = 'PF_condom_'
  return (
    # across years (same type)
    P[k+'msp_2006'] < P[k+'msp_2016'] and
    P[k+'cas_2006'] < P[k+'cas_2016'] and
    P[k+'swo_2002'] < P[k+'swo_2011'] < P[k+'swo_2014'] and
    P[k+'swr_2002'] < P[k+'swr_2011'] < P[k+'swr_2014'] and
    # across types (same year)
    P[k+'msp_2006'] < P[k+'cas_2006'] and
    P[k+'msp_2016'] < P[k+'cas_2016'] and
    P[k+'swr_2002'] < P[k+'swo_2002'] and 
    P[k+'swr_2011'] < P[k+'swo_2011'] and 
    P[k+'swr_2014'] < P[k+'swo_2014']
  )

def get_circumcision(P): # [OK]
  # Rbeta_circum.shape = (a:2, s:2)
  # (SHIMS2), SDHS2006, Bicego2013, MICS2014, SHIMS2, "COP20", assume
  PF_circum_t = ta.tarray(
     [1980.0,2006.5,2011.0,2014.8,2016.5,2020.0,2050,2051],
     [  .007,  .082,  .171,  .250,  .300,  .370,*2*[P['PF_circum_2050']]]
  ).reshape([1,1,1,1,1,1,1,1])
  Rbeta_circum =  np.array([ # [women,men]
      [1,.50], # vaginal - Boily2009,Hughes2012,Patel2014
      [1,.27], # anal - Wiysonge2011
    ]).reshape([2,1,2,1,1,1,1,1])
  return {
    'PF_circum_t': PF_circum_t,
    'Rbeta_circum': Rbeta_circum,
  }

def get_mix(P): # [OK]
  # pref_pii.shape = (p:4, i:4, i':4)
  pref_pii = np.zeros((4,4,4))
  pref_pii[0,0,0] = P['pref_msp_xl']
  pref_pii[0:2,2:,2:] = P['pref_mcx_swx']
  return {
    'pref_pii': pref_pii,
    'mix': np.zeros((4,2,4,2,4)), # initialize
    'mix_mask': np.ones((4,2,4,2,4)), # TODO
    't0_tpaf': np.inf,
  }

def get_mix_mask(mask=None,p=None,sfr=None,ifr=None,sto=None,ito=None):
  if mask is None: mask = np.ones((4,2,4,2,4))
  if p   is None: p   = slice(None)
  if sfr is None: sfr = slice(None)
  if ifr is None: ifr = slice(None)
  if sto is None: sto = slice(None)
  if ito is None: ito = slice(None)
  mask[p,sfr,ifr,sto,ito] = 0
  return mask

# HIV ----------------------------------------------------------------------------------------------

def get_hiv_prog(P): # [OK]
  # h: sus, acute, >500, <500, <350, <200 (AIDS)
  # c: undiag, diag, unlinked, on art, vls
  dur_h  = np.array([P['dur_acute'],3.5-P['dur_acute'],3.74,5.26]).reshape([1,1,1,4,1]) # Lodi2011,Mangal2017
  prog_h = 1/dur_h
  unprog_h = np.array([[2,.1],[2,.1],[1.5,.1]]).reshape([1,1,1,3,2]) # Battegay2006,Lawn2006
  death_h = np.array([0,0,.004,.02,.04,.20]).reshape([1,1,1,6,1]) # Badri2006,Anglaret2012,Mangal2017
  Rdeath_c = np.array([1,1,1,.25,.5]).reshape([1,1,1,1,5]) # Gabillard2013,Lundgren2015
  return {
    'dur_h': dur_h,
    'prog_h': prog_h,
    'unprog_h': unprog_h,
    'death_hc': death_h * Rdeath_c,
  }

# cascade ------------------------------------------------------------------------------------------

def get_diag(P):
  t_dx = [1980,1990,2002,2006,2011,2016,2051]
  dx_wq_t = np.array([0,0,P['dx_w_2002'],P['dx_w_2006'],P['dx_wq_2011'],
    *2*[P['dx_wq_2011']*(1+P['aRdx_wq_16:11'])]])
  Rdx_wq_t  = np.array([1,1,1,1,1,1,1]) # dummy
  Rdx_fsw_t = np.array([1,1,1,1,1+P['aRdx_fsw:wq_2011'],*2*[1+P['aRdx_fsw:wq_2016']]])
  Rdx_m_t   = np.array([1,1,0.1,P['Rdx_m:w_2006'],*3*[P['Rdx_m:wq_2011']]])
  dx_sit = ta.tarray(t_dx,P['Rdx_global']*dx_wq_t*np.array(
    [[Rdx_wq_t,Rdx_wq_t,Rdx_fsw_t,Rdx_fsw_t],
     [ Rdx_m_t, Rdx_m_t,  Rdx_m_t,  Rdx_m_t]])).reshape((2,4,1,1))
  return {
    'dx_sit': dx_sit,
  }

def get_treat(P):
  Rtx_si = np.array([[1,1,P['Rtx_fsw:wq'],P['Rtx_fsw:wq']],[1,1,1,1]]).reshape([2,4,1,1,1])
  tx_sit = ta.tarray([1980,2003,2010,2012,2018,2051],
    Rtx_si * np.array([0,0,P['tx_2010'],P['tx_2012'],12,12]).reshape([1,1,1,1,6]))
  Rtx_ht = ta.tarray(
     [1980,2003,2004,2010,2011,2015,2016,2017,2018,2051], [
     [   0,   0,   0,   0,   0,   0,   0,   0,   1,   1], # acute:           scale-up 2017-2018
     [   0,   0,   0, .05, .05, .05, .05, .15,   1,   1], # cd4 > 500:       scale-up 2017-2018
     [   0,   0, .15, .15, .15, .15,   1,   1,   1,   1], # 350 < cd4 < 500: scale-up 2015-2016
     [   0,   0, .35, .35,   1,   1,   1,   1,   1,   1], # 200 < cd4 < 350: scale-up 2010-2011
     [   0,   0,   1,   1,   1,   1,   1,   1,   1,   1], # cd4 < 200:       scale-up 2003-2004
  ]).reshape([1,1,1,5])
  vx = np.array(1/P['ivx'])
  unvx_t = ta.tarray([1980,2010,2018,2051],[.15,.15,.05,.05])
  Runvx_si = np.array([
    [1,1,P['Runvx_fsw:wq'],P['Runvx_fsw:wq']],
    [P['Runvx_m:wq'],P['Runvx_m:wq'],P['Runvx_m:wq'],P['Runvx_m:wq']]]).reshape((2,4,1,1))
  revx_t = ta.tarray([1980,2010,2018,2051],P['revx_2010']*np.array([1,1,1.5,1.5]))
  return {
    'tx_sit': tx_sit,
    'Rtx_ht': Rtx_ht,
    'vx': vx,
    'unvx_t': unvx_t,
    'Runvx_si': Runvx_si,
    'revx_t': revx_t,
  }

def get_scen(P):
  return {
    'Rdx_scen': np.ones((2,4,1,1)),
    'Rtx_scen': np.ones((2,4,1,1)),
    'Rux_scen': np.ones((2,4,1,1)),
  }
