# functions to specify & sample model parameters (inputs)

import numpy as np
from copy import copy
from scipy.optimize import nnls
from utils import _,NAN,log,stats,flatten,dict_split,linear_comb,interval_qs
from utils import tarray as ta
from model import tol

# ------------------------------------------------------------------------------
# sampling functions

def def_sample_distrs():
  # define a dict (D) of param sampling distribs
  return {
  # popuation
  't0_hiv':               stats.uniform(l=1980,h=1985),
  'PX_w_fsw':             stats.betabin(p=.0288,n=121),
  'PX_w_h':               stats.betabin(p=.1778,n=66),
  'PX_m_m':               stats.betabin(p=.1331,n=360),
  'dur_fsw':              stats.gamma(m=4.065,sd=1.039),
  'dur_sw_h':             stats.gamma(m=0.500,sd=0.216),
  'dur_cli':              stats.gamma(m=8.626,sd=2.832),
  'turn_xm_xl':           stats.gamma(m=.2156,sd=.1177),
  'Pturn_fsw_m:l':        stats.betabin(p=.724,n=18),
  'Pturn_cli_m:l':        stats.betabin(p=.602,n=7),
  'growth_2050':          stats.uniform(l=.007,h=.015),
  # partnership rates
  'C12m_msp_xl':          stats.betabin(p=.370,n=57),
  'C12m_cas_xl':          stats.betabin(p=.366,n=28),
  'C12m_cas_wm':          stats.gamma(m=1.575,sd=0.204),
  'RC_cas_cli:wm':        stats.uniform(l=.25,h=1),
  'C1m_swo_fsw_l':        stats.gamma(m=3.5,sd=0.7),
  'C1m_swr_fsw_l':        stats.gamma(m=6.0,sd=1.2),
  'C1m_swo_fsw_h':        stats.gamma(m=14,sd=2.8),
  'C1m_swr_fsw_h':        stats.gamma(m=21,sd=4.2),
  'KF_swx_cli':           stats.gamma(m=40.5,sd=13.9),
  'RKF_swx_cli_h:l':      stats.gamma(m=2.03,sd=0.230),
  # sex / mixing
  'F_msp':                stats.gamma(m=77.3,sd=33.7),
  'RF_cas:msp':           stats.uniform(l=.5,h=2),
  'dur_msp':              stats.uniform(l=14.5,h=18.5),
  'dur_cas':              stats.gamma(m=.743,sd=0.324),
  'dur_swr':              stats.gamma(m=0.495,sd=0.216),
  'F_swr':                stats.uniform(l=12,h=36),
  'PF_ai_mcx':            stats.gamma(m=.0573,sd=.0423),
  'RPF_ai_swx:mcx':       stats.uniform(l=1,h=2),
  'lpref_msp_xl':         stats.uniform(l=-2,h=+2),
  'lpref_mcx_swx':        stats.uniform(l=-2,h=+2),
  # condom use
  'Rbeta_condom':         stats.betabin(p=.266,n= 32),
  'RPF_condom_a:v':       stats.betabin(p=.768,n= 12),
  'RPF_condom_1996':      stats.uniform(l=0,h=1),
  'PF_condom_msp_2006':   stats.betabin(p=.230,n=100),
  'PF_condom_msp_2016':   stats.betabin(p=.416,n= 75),
  'PF_condom_cas_2006':   stats.betabin(p=.598,n=100),
  'PF_condom_cas_2016':   stats.betabin(p=.694,n=100),
  'PF_condom_swo_2002':   stats.betabin(p=.432,n=  9),
  'PF_condom_swo_2011':   stats.betabin(p=.777,n= 21),
  'PF_condom_swo_2014':   stats.betabin(p=.787,n= 14),
  'PF_condom_swr_2002':   stats.betabin(p=.337,n= 13),
  'PF_condom_swr_2011':   stats.betabin(p=.754,n= 24),
  'PF_condom_swr_2014':   stats.betabin(p=.759,n= 11),
  'PF_circum_2050':       stats.betabin(p=.724,n= 18),
  # transmit prob
  'beta_0':               stats.gamma(m=.00131,sd=.00052),
  'Rbeta_acute':          stats.gamma(m=6.015,sd=3.653),
  'Rbeta_350':            stats.gamma(m=1.586,sd=0.1532),
  'Rbeta_200':            stats.gamma(m=8.203,sd=2.182),
  'Rbeta_vi_rec':         stats.uniform(l=1,h=2),
  'aRbeta_gud_sus':       stats.gamma(m=2.05,sd=1.546),
  'aRbeta_gud_inf':       stats.gamma(m=0.99,sd=0.577),
  'dur_acute':            stats.gamma(m=.172,sd=.1285),
  'P_gud_fsw_l':          stats.betabin(p=.295,n=78),
  'RP_gud_fsw_h:l':       stats.gamma(m=1.283,sd=0.153),
  'RP_gud_2030':          stats.uniform(l=.2,h=1),
  'iP_gud_h:l':           stats.uniform(l=0,h=1),
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

def get_lp(D,P):
  # log prior probability of observing P under D, including constraints
  return sum(D[k].logpdf(P[k]) for k in D.keys()) + \
    sum(0 if constr(P) else -np.inf for constr in def_constrs())

def get_n_sample_lhs(D,n,seed):
  # get n samples (list of dicts) via Latin Hypercube Sampling (LHS)
  # from the distribs D (dict)
  Qs = stats.lhs(len(D),n,seed) # Qs is (n x len(D)) array of quantiles
  return [{k:D[k].ppf(qk) for k,qk in zip(D,Q)} for Q in Qs]

def get_n_sample_constr(D,seeds,Ps=None):
  # get n samples (list of dicts) via get_random_constr
  if Ps is None: Ps = [{} for seed in seeds]
  return [get_sample_constr(D,seed,P) for seed,P in zip(seeds,Ps)]

def get_sample_constr(D,seed=None,P=None):
  # get param set (dict) via random sampling with constraints
  if seed is not None: np.random.seed(seed)
  if P is None: P = {}
  P.update({k:Dk.rvs() for k,Dk in D.items()}) # initial samples
  P['id'] = seed
  for constr,ks in def_constrs().items():
    resample_constr(D,P,constr,ks)
  return P

def resample_constr(D,P,constr,ks):
  # resample the list of params (ks) randomly until constr(P) = True
  while not constr(P):
    for k in ks:
      P[k] = D[k].rvs()
  return P

def def_constrs():
  # define constraints (dict) like {fun: [keys], ...} for resample_constr
  k = 'PF_condom_' # convenience
  return {
    constr_F:      ['F_swr'],
    constr_acute:  ['Rbeta_acute','dur_acute'],
    constr_gud:    ['P_gud_fsw_l','RP_gud_fsw_h:l'],
    constr_condom: [k+'msp_2006',k+'msp_2016',
                    k+'cas_2006',k+'cas_2016',
                    k+'swo_2002',k+'swo_2011',k+'swo_2014',
                    k+'swr_2002',k+'swr_2011',k+'swr_2014'],
  }

def get_depend(P,**kwds):
  # add dependent params to dict P, which already contains sampled params
  P.setdefault('foi_mode','base')
  P.update(kwds)
  # order below matters for some dependencies
  P.update(get_F(P))
  P.update(get_PX(P))
  P.update(get_birth_death(P))
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

def get_n_all(seeds,**kwds):
  # get a list (Ps) of n complete param dicts (P) for given seeds
  log(2,'params.get_n_all: N = '+str(len(seeds)))
  D   = def_sample_distrs()
  Dc  = dict_split(D,flatten(def_constrs().values())) # split un/constrained
  Ps0 = get_n_sample_lhs(D,len(seeds),seed=seeds[0])  # unconstrained
  Ps  = get_n_sample_constr(Dc,seeds,Ps=Ps0)          # constrained
  return [get_depend(P,**kwds) for P in Ps]

# ------------------------------------------------------------------------------
# population

def get_PX(P):
  # PX: risk group sizes (relative); N = (absolute, '000s)
  # PX.shape = (s:2, i:4, k:5, h:6, c:5)
  NX0   = np.array(243) # REF: WorldBank for Eswatini 1980 pop ('000s)
  PX0_k = np.array([1,0,0,0,0]).reshape([1,1,5,1,1])   # initial EPA
  PX0_h = np.array([1,0,0,0,0,0]).reshape([1,1,1,6,1]) # initial health
  PX0_c  = np.array([1,0,0,0,0]).reshape([1,1,1,1,5])  # initial care
  PX_h_hiv = np.array([0,5,65,30,0,0])*1e-6 # seed HIV; REF: assume
  PX_h_hiv[0] = 1 - PX_h_hiv.sum() # susceptible
  P = get_PX_swx(P)
  PX_si = np.zeros((2,4)) # risk (sex/activity) groups
  # FSW
  PX_si[0,2] = P['PX_fsw'] * (1-P['PX_fsw_h'])
  PX_si[0,3] = P['PX_fsw'] * P['PX_fsw_h']
  # clients
  PX_si[1,3] = P['PX_cli'] * P['PX_cli_h']
  PX_si[1,2] = P['PX_cli'] * (1-P['PX_cli_h'])
  # non-client men
  PX_si[1,1] = P['PX_m_m'] * (1-P['PX_w'])
  PX_si[1,0] = 1 - P['PX_w'] - PX_si[1,1:].sum()
  # non-FSW women
  PX_si[0,1] = P['PX_w'] * P['PX_w_h'] - P['PX_fsw']
  PX_si[0,0] = P['PX_w'] - PX_si[0,1:].sum()
  return {
    'PX_m_cli': P['PX_cli'] / (1 - P['PX_w']), # % clients / men
    'PX_s':  PX_si.sum(axis=1),                # % sex / overall
    'PX_si': PX_si,                            # % activity & sex / pop
    'PX_si_s': PX_si / PX_si.sum(axis=1)[:,_], # % activity / sex
    'X0': NX0 * PX_si[:,:,_,_,_] * PX0_k * PX0_h * PX0_c, # full pop (t=0)
    'PX_h_hiv': PX_h_hiv, # seed HIV
  }

def get_PX_swx(P):
  P['PX_w']     = .52 # % women / overall; REF: WorldBank
  P['PX_fsw_h'] = .2  # % high risk / FSW; REF: assume
  P['PX_cli_h'] = .2  # % high risk / cli; REF: assume
  P['F_swo']    = 12  # see get_F()
  # calculating client size from total sex work acts
  P['PX_fsw']   = P['PX_w'] * P['PX_w_fsw']
  P['XKF_swo']  = P['PX_fsw'] * P['C2K_p'][2] * P['F_swo'] * \
    linear_comb(P['PX_fsw_h'],P['C1m_swo_fsw_h'],P['C1m_swo_fsw_l'])
  P['XKF_swr']  = P['PX_fsw'] * P['C2K_p'][2] * P['F_swr'] * \
    linear_comb(P['PX_fsw_h'],P['C1m_swr_fsw_h'],P['C1m_swr_fsw_l'])
  # XKF = total sex acts (pop x partners x sex freq)
  P['PX_cli'] = (P['XKF_swo'] + P['XKF_swr']) / P['KF_swx_cli']
  # note: KF_swx_cli is true client average: see get_K()
  return P

def get_birth_death(P):
  # birth rate: age 15+ entry; death rate: non-HIV mortality
  death = np.array(1/35 + (1-.64)*.0144) # death rate ~ .034 for 15-49 years + non-HIV mortality
  birth = ta.tarray([1980,2000,2010,2020,2050],np.array([.04,.03,.015,.015,P['growth_2050']])+death)
  return {
    'birth_t': birth,
    'death': death,
  }

def get_turnover(P):
  # define A,b for: A*w = b (linear constraints), where "w" = [e,t]:
  # ei: size of group i among entrants; tij: turnover from i to j
  #  0, 1, 2, 3,   4,  5,  6,   7,  8,  9,  10, 11, 12,  13, 14, 15
  # e1,e2,e3,e4, t01,t02,t03, t10,t12,t13, t20,t21,t23, t30,t31,t32
  # can't solve turnover here, because birth_t changes each timestep
  PXe  = np.zeros((2,4))      # model entrants "e" (s,i)
  turn = np.zeros((2,4,4))    # turnover matrix "t" (s,i.from,i.to)
  b4   = np.eye(4,dtype=bool) # boolean I matrix (pre-compute for speed)
  A = np.zeros((2,16,16))     # constraint matrix "A"
  b = np.zeros((2,16))        # constraint vector "b"
  for s in (0,1): # note: {variables} = NAN below computed at runtime in solve_turnover
    x = P['PX_si'][s,:] # steady-state group size
    p = (P['PX_fsw_h'],P['PX_cli_h'])[s] # 0.20 (assumed)
    dur = (P['dur_fsw'],P['dur_cli'])[s] # duration for FSW or clients
    # steady-state constraints
    b[s,0:4] = NAN # x*{birth_t}
    A[s,0,:] = (NAN,0,0,0,-x[0],-x[0],-x[0],+x[1],    0,    0,+x[2],    0,    0,+x[3],    0,    0)
    A[s,1,:] = (0,NAN,0,0,+x[0],    0,    0,-x[1],-x[1],-x[1],    0,+x[2],    0,    0,+x[3],    0)
    A[s,2,:] = (0,0,NAN,0,    0,+x[0],    0,    0,+x[1],    0,-x[2],-x[2],-x[2],    0,    0,+x[3])
    A[s,3,:] = (0,0,0,NAN,    0,    0,+x[0],    0,    0,+x[1],    0,    0,+x[2],-x[3],-x[3],-x[3])
    # other constraints
    A[s,4,0:4],b[s,4] = 1, x.sum() # e.sum = x.sum
    A[s,5,1],b[s,5] = 1, x[1]*1.0 # e1 = x1
    A[s,6,2],b[s,6] = 1, NAN # e2 = (x2+x3)*{R}
    A[s,7,3],b[s,7] = 1, 0.0 # e3 = 0
    A[s,8,10:12],b[s,8] = 1, (1/dur - P['death'])/(1-p) # dur 2 (adj for +3)
    A[s,9,15],   b[s,9] = 1, (1/P['dur_sw_h'] - P['death']) # dur 3
    A[s,10, 7],b[s,10] = 1, P['turn_xm_xl'] # t10
    A[s,11, 6],b[s,11] = 1, 0.0 # t03 = 0
    A[s,12, 9],b[s,12] = 1, 0.0 # t13 = 0
    A[s,13,13],b[s,13] = 1, 0.0 # t30 = 0
    A[s,14,14],b[s,14] = 1, 0.0 # t31 = 0
    A[s,15, :],b[s,15] = (0,0,0,0, 0,0,0, 0,0,0, 0,0,1-p, 0,0,-p), p # t23*x2 - t32*x3 = x3*{birth_t}
  # print(np.array(sympy.Matrix(A[s,:,:]).rref()[0],dtype=float).round(3)) # DEBUG @ NAN = 0.03
  # print(np.linalg.matrix_rank(A)) # DEBUG ~= [16,16] @ NAN = 0.03
  return {
    'PXe_si': PXe,
    'turn_sii': turn,
    'turn_A_s': A,
    'turn_b_s': b,
    'b4': b4,
  }

def solve_turnover(P,t):
  # replace NAN with runtime values (depend on birth_t) & solve "w" in A*w = b
  v = P['birth_t'](t)
  P['turn_b_s'][:,0:4]     = v * P['PX_si']
  P['turn_A_s'][:,0:4,0:4] = v * P['b4']
  for s,R,dur in zip((0,1),(2.0,1.5),(P['dur_fsw'],P['dur_cli'])):
    P['turn_b_s'][s,6]  = P['PX_si'][s,2:].sum() * np.minimum(R, (v - P['death'] + 1/dur) / v)
    P['turn_b_s'][s,15] = v * P['PX_fsw_h']
    w,err = nnls(P['turn_A_s'][s],P['turn_b_s'][s])  # non-negative least squares
    P['PXe_si'][s] = w[0:4]                          # e: entrant group sizes
    P['turn_sii'][s,np.logical_not(P['b4'])] = w[4:] # t: turnover matrix
    if err > tol:
      P['PXe_si'][s] = -np.inf # then system.solve will fail (gracefully)
      # raise Exception('Cannot solve turnover (s={},t={},seed={}) error: {}'.format(s,t,P['id'],err))
  return v,P['PXe_si'],P['turn_sii']

# ------------------------------------------------------------------------------
# force of infection

def get_beta_a(P):
  # beta: probability of transmission per sex act
  # beta_a.shape = (a:2, p:4, s:2, i:4, s':2, i':4, h':6, c':5, (t))
  Rbeta_ar = 10 # ar = anal receptive
  Rbeta_as = np.array([[P['Rbeta_vi_rec'],1],[Rbeta_ar,1]]).reshape([2,1,2,1,1,1,1,1])
  Rbeta_as = Rbeta_as / Rbeta_as[0,:].mean() # as = act, sex
  # prevalence of GUD (genital ulcer disease)
  P_gud_0     = .07 # REF: SDHS2006
  P_gud_m     = linear_comb(P['iP_gud_h:l'],P_gud_0,P['P_gud_fsw_l'])
  P_gud_cli_l = linear_comb(.5,P_gud_m,P['P_gud_fsw_l'])
  P_gud_cli_h = P['P_gud_fsw_l']
  P_gud = np.array([ # TODO: rename P_gud_si (?)
    [P_gud_0, P_gud_m, P['P_gud_fsw_l'], P['P_gud_fsw_l']*P['RP_gud_fsw_h:l']], # MAN
    [P_gud_0, P_gud_m, P_gud_cli_l, P_gud_cli_h], # MAN
  ])
  RP_gud_t = ta.tarray([1980,2000,2010,2030,2051],[1,1,1,*2*[P['RP_gud_2030']]])
  # relative beta: health & care
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

def constr_acute(P):
  return 1 <= P['Rbeta_acute'] * P['dur_acute'] <= 63 # REF: Bellan2015

def constr_gud(P):
  return (
    P['P_gud_fsw_l'] > .07 and # REF: SDHS2006
    P['P_gud_fsw_l'] * P['RP_gud_fsw_h:l'] < 1
  )

def get_F(P):
  # F: frequency of sex per-partnership, dur: partnership duration
  # F_ap.shape = (a:2, p:4)
  PF_ai = P['PF_ai_mcx'] * np.array([1,1,P['RPF_ai_swx:mcx'],P['RPF_ai_swx:mcx']])
  F_p   = np.array([ P['F_msp'], P['F_msp']*P['RF_cas:msp'], 12, P['F_swr'] ])
  F_ap  = F_p.reshape([1,4]) * np.array([1-PF_ai,PF_ai]).reshape([2,4])
  dur_p = np.array([ P['dur_msp'], P['dur_cas'], 1/12, P['dur_swr'] ])
  dur_p_1 = np.minimum(dur_p,1) # for foi.cases = ('ry','py')
  # C2K = factor to convert reported (C) -> current (K) given recall periods (in years)
  C2K_p = dur_p / (dur_p + np.array([1, 1, 1/12, 1/12]))
  return {
    'F_ap': F_ap,
    'dur_p': dur_p,
    'dur_p_1': dur_p_1,
    'C2K_p': C2K_p,
  }

def constr_F(P):
  C2K_swo = 1/2 # (1/12) / (1/12 + 1/12)
  C2K_swr = P['dur_swr'] / (P['dur_swr'] + 1/12)
  return (
    P['C1m_swo_fsw_h'] * C2K_swo * 12 + \
    P['C1m_swr_fsw_h'] * C2K_swr * P['F_swr']
  ) < 2*365 # max sex work partners: 2 per day; REF: assume

def get_K(P):
  # K: current partner numbers
  # K.shape = (p:4, s:2, i:4)
  PX_si = P['PX_si']
  F = np.squeeze(P['F_ap'].sum(axis=0))
  K_psi = np.zeros((4,2,4,1)) # TODO: need singleton?
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
  K_psi[2,0,3] = P['C2K_p'][2] * P['C1m_swo_fsw_h']
  K_psi[3,0,2] = P['C2K_p'][3] * P['C1m_swr_fsw_l']
  K_psi[3,0,3] = P['C2K_p'][3] * P['C1m_swr_fsw_h']
  # swx: clients
  wPX = (PX_si[1,2] + P['RKF_swx_cli_h:l'] * PX_si[1,3])
  K_psi[2,1,2] = P['XKF_swo'] / F[2] / wPX
  K_psi[2,1,3] = P['XKF_swo'] / F[2] / wPX * P['RKF_swx_cli_h:l']
  K_psi[3,1,2] = P['XKF_swr'] / F[3] / wPX
  K_psi[3,1,3] = P['XKF_swr'] / F[3] / wPX * P['RKF_swx_cli_h:l']
  KF_psi = K_psi.sum(axis=3) * F[:,_,_]
  # aK_pk: K reduction for ptr-type "p" in EPA stratum "k"
  aK_pk = np.array([[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]]).reshape([4,1,1,5])
  aK_pk = aK_pk * (K_psi > 0) # don't reduce partners that don't exist
  return {
    'K_psi': K_psi,
    'KF_psi': KF_psi,
    'aK_pk': aK_pk,
  }

def get_condom(P):
  k = 'PF_condom_' # convenience
  R = P['RPF_condom_1996'] # allow more variation vs linear interp to 2002/06
  PF_condom_t = ta.tarray([1980,1988,1996,2002,2006,2011,2014,2016,2050],
    [[0,  0,R*P[k+'msp_2006'],NAN,P[k+'msp_2006'],NAN,NAN,*2*[P[k+'msp_2016']] ], # main
     [0,  0,R*P[k+'cas_2006'],NAN,P[k+'cas_2006'],NAN,NAN,*2*[P[k+'cas_2016']] ], # casual
     [0,NAN,R*P[k+'swo_2002'],P[k+'swo_2002'],NAN,P[k+'swo_2011'],*3*[P[k+'swo_2014']] ], # sw-new
     [0,NAN,R*P[k+'swr_2002'],P[k+'swr_2002'],NAN,P[k+'swr_2011'],*3*[P[k+'swr_2014']] ]] # sw-reg
  ).reshape([1,4,1,1,1,1,1,1])
  RPF_condom_a = np.array([1,P['RPF_condom_a:v']]).reshape([2,1,1,1,1,1,1,1])
  return {
    'PF_condom_t': PF_condom_t,
    'RPF_condom_a': RPF_condom_a,
  }

def constr_condom(P):
  # TODO: (?) NERCHA2018: cas condom use may be declining
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

def get_circumcision(P):
  # Rbeta_circum.shape = (a:2, s:2)
  # REF: (SHIMS2), SDHS2006, Bicego2013, MICS2014, SHIMS2, "COP20", assume
  PF_circum_t = ta.tarray(
     [1980.0,2006.5,2011.0,2014.8,2016.5,2020.0,2050,2051],
     [  .007,  .082,  .171,  .250,  .300,  .370,*2*[P['PF_circum_2050']]]
  ).reshape([1,1,1,1,1,1,1,1])
  Rbeta_circum =  np.array([ # women, men
      [1,.50], # vaginal; REF: Boily2009,Hughes2012,Patel2014
      [1,.27], # anal; REF Wiysonge2011
    ]).reshape([2,1,2,1,1,1,1,1])
  return {
    'PF_circum_t': PF_circum_t,
    'Rbeta_circum': Rbeta_circum,
  }

def get_mix(P):
  # pref_pii: log-likelihood of mixing of i with i' for ptr-type p
  # pref_pii.shape = (p:4, i:4, i':4)
  pref_pii = np.zeros((4,4,4))
  pref_pii[0,0,0] = np.exp(P['lpref_msp_xl'])
  pref_pii[0:2,2:,2:] = np.exp(P['lpref_mcx_swx'])
  return {
    'pref_pii': pref_pii,
    'mix': np.zeros((4,2,4,2,4)), # initialize
    'mix_mask': np.ones((4,2,4,2,4)), # for tpaf
    't0_tpaf': np.inf,
  }

def get_mix_mask(mask=None,p=None,fs=None,fi=None,ts=None,ti=None):
  # disable (mask) transmission via selected pathways;
  # p: ptr-type; f*: from; t*: to; *s: sex; *i: activity
  if mask is None: mask = np.ones((4,2,4,2,4))
  if p  is None: p  = slice(None)
  if fs is None: fs = slice(None)
  if fi is None: fi = slice(None)
  if ts is None: ts = slice(None)
  if ti is None: ti = slice(None)
  mask[p,ts,ti,fs,fi] = 0
  return mask

# ------------------------------------------------------------------------------
# HIV

def get_hiv_prog(P):
  # h: sus, acute, >500, <500, <350, <200 (AIDS)
  # c: undiag, diag, unlinked, on art, vls
  dur_h  = np.array([P['dur_acute'],3.5-P['dur_acute'],3.74,5.26]).reshape([1,1,1,4,1]) # REF: Lodi2011,Mangal2017
  prog_h = 1/dur_h # progression rate = 1 / mean duration
  unprog_h = np.array([[2,.1],[2,.1],[1.5,.1]]).reshape([1,1,1,3,2]) # REF: Battegay2006,Lawn2006
  death_h = np.array([0,0,.004,.02,.04,.20]).reshape([1,1,1,6,1]) # REF: Badri2006,Anglaret2012,Mangal2017
  Rdeath_c = np.array([1,1,1,.25,.5]).reshape([1,1,1,1,5]) # REF: Gabillard2013,Lundgren2015
  return {
    'dur_h': dur_h,
    'prog_h': prog_h, # untreated HIV progression
    'unprog_h': unprog_h, # treated HIV
    'death_hc': death_h * Rdeath_c, # HIV death
  }

# ------------------------------------------------------------------------------
# cascade

def get_diag(P):
  # dx: diagnosis rates (by group & time)
  t_dx = [1980,1990,2002,2006,2011,2016,2051]
  dx_wq_t = np.array([0,0,P['dx_w_2002'],P['dx_w_2006'],P['dx_wq_2011'],
    *2*[P['dx_wq_2011']*(1+P['aRdx_wq_16:11'])]])
  Rdx_wq_t  = np.array([1,1,1,1,1,1,1]) # dummy
  Rdx_fsw_t = np.array([1,1,1,1,1+P['aRdx_fsw:wq_2011'],*2*[1+P['aRdx_fsw:wq_2016']]])
  Rdx_m_t   = np.array([1,1,0.1,P['Rdx_m:w_2006'],*3*[P['Rdx_m:wq_2011']]])
  dx_sit = ta.tarray(t_dx,P['Rdx_global']*dx_wq_t*np.array(
    [[Rdx_wq_t,Rdx_wq_t,Rdx_fsw_t,Rdx_fsw_t],
     [ Rdx_m_t, Rdx_m_t,  Rdx_m_t,  Rdx_m_t]])).reshape([2,4,1,1])
  return {
    'dx_sit': dx_sit,
  }

def get_treat(P):
  # tx: treatment rates (by group, HIV stage, & time)
  Rtx_si = np.array([[1,1,P['Rtx_fsw:wq'],P['Rtx_fsw:wq']],[1,1,1,1]]).reshape([2,4,1,1,1])
  tx_sit = ta.tarray([1980,2003,2010,2012,2018,2051],
    Rtx_si * np.array([0,0,P['tx_2010'],P['tx_2012'],12,12]).reshape([1,1,1,1,6]))
  Rtx_ht = ta.tarray(
     [1980,2003,2004,2010,2011,2015,2016,2017,2018,2051],
    [[   0,   0,   0,   0,   0,   0,   0,   0,   1,   1], # acute:           scale-up 2017-2018
     [   0,   0,   0, .05, .05, .05, .05, .15,   1,   1], # cd4 > 500:       scale-up 2017-2018
     [   0,   0, .15, .15, .15, .15,   1,   1,   1,   1], # 350 < cd4 < 500: scale-up 2015-2016
     [   0,   0, .35, .35,   1,   1,   1,   1,   1,   1], # 200 < cd4 < 350: scale-up 2010-2011
     [   0,   0,   1,   1,   1,   1,   1,   1,   1,   1]] # cd4 < 200:       scale-up 2003-2004
  ).reshape([1,1,1,5])
  # vx: viral suppression; unvx: treatment fail/discontinue; revx: re-suppress
  vx = np.array(1/P['ivx'])
  Runvx_si = np.array([
    [1,1,P['Runvx_fsw:wq'],P['Runvx_fsw:wq']],
    [P['Runvx_m:wq'],P['Runvx_m:wq'],P['Runvx_m:wq'],P['Runvx_m:wq']]]).reshape([2,4,1,1,1])
  unvx_sit = ta.tarray([1980,2010,2018,2051],
    Runvx_si * np.array([.15,.15,.05,.05]).reshape([1,1,1,1,4]))
  revx_t = ta.tarray([1980,2010,2018,2051],P['revx_2010']*np.array([1,1,1.5,1.5]))
  return {
    'tx_sit': tx_sit,
    'Rtx_ht': Rtx_ht,
    'vx': vx,
    'unvx_sit': unvx_sit,
    'revx_t': revx_t,
  }

def get_scen(P):
  # relative rates (dx, tx, unvx) for scenarios
  return {
    'Rdx_scen': np.ones((2,4,1,1)),
    'Rtx_scen': np.ones((2,4,1,1)),
    'Rux_scen': np.ones((2,4,1,1)),
  }
