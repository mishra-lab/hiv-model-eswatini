import numpy as np
from copy import copy
from utils import _,NAN,stats,flatten,dict_split,linear_comb,interval_qs,log
from utils import tarray as ta

# main -------------------------------------------------------------------------

def get_all(P,seed=None,**kwds):
  P = P if P is not None else \
      get_sample_random(seed=seed)
  P.update(kwds)
  # order matters for some dependencies
  P.update(get_F(P))
  P.update(get_PX(P))
  P.update(get_birth_death(P))
  P.update(get_turnover(P))
  P.update(get_C(P))
  P.update(get_mix(P))
  P.update(get_beta_a(P))
  P.update(get_condom(P))
  P.update(get_circumcision(P))
  P.update(get_hiv_prog(P))
  P.update(get_diag(P))
  P.update(get_treat(P))
  # check.ch_pn(P) # TEMP
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
    check_F:      ['PF_ai_mcq','PF_ai_swq'],
    check_acute:  ['Rbeta_acute','dur_acute'],
    check_gud:    ['P_gud_fsw_l','RP_gud_fsw_h:l'],
    check_condom: [k+'msp_1988',k+'msp_2006',k+'msp_2016',k+'cas_1988',k+'cas_2006',k+'cas_2016',
                   k+'swo_2002',k+'swo_2011',k+'swo_2014',k+'swr_2002',k+'swr_2011',k+'swr_2014']
  }

def def_sample_distrs():
  return {
  # PX
  't0_hiv':               stats.uniform(l=1980,h=1985),
  'PX_w_fsw':             stats.beta_binom(p=.028,n=112),
  'PX_w_h':               stats.beta_binom(p=.177,n=75),
  'PX_m_m':               stats.beta_binom(p=.133,n=133),
  'dur_fsw_l':            stats.gamma_p(p=3.6,v=.96),
  'dur_fsw_h':            stats.gamma_p(p=10,v=0.26),
  'dur_cli':              stats.gamma_p(p=10,v=5.32),
  'Pturn_fsw_h:m':        stats.beta_binom(p=.724,n=18),
  'LORturn_wq_sus:hiv':   stats.gamma_p(p=1,v=.605),
  'growth_2050':          stats.uniform(l=.007,h=.015),
  # C
  'C_swo_fsw_l':          stats.gamma_p(p=49.1,v=115),
  'C_swr_fsw_l':          stats.gamma_p(p=5.7,v=3.32),
  'RC_swo_fsw_h:l':       stats.gamma_p(p=2.0,v=.053),
  'RC_swr_fsw_h:l':       stats.gamma_p(p=1.5,v=.01),
  'CF_swq_cli':           stats.gamma_p(p=60,v=197),
  'RCF_swq_cli_h:l':      stats.gamma_p(p=2.0,v=.05),
  'C_cas_am':             stats.gamma_p(p=2.0,v=.213),
  'RC_cas_fsw:am':        stats.uniform(l=.25,h=1),
  'RC_cas_cli:am':        stats.uniform(l=.25,h=1),
  # F
  'F_mcq':                stats.gamma_p(p=78,v=1130),
  'dur_msp':              stats.gamma_p(p=18,v=11),
  'dur_cas':              stats.gamma_p(p=.25,v=1.16e-2),
  'F_swr':                stats.uniform(l=12,h=48),
  'PF_ai_mcq':            stats.beta_binom(p=.059,n=30),
  'PF_ai_swq':            stats.beta_binom(p=.097,n=14),
  # mixing
  'pref_msp_al':          stats.gamma_p(p=2,v=.2),
  'pref_msp_asw':         stats.gamma_p(p=3,v=.8),
  # condoms
  'RPF_condom_a:v':       stats.beta_binom(p=.768,n= 12),
  'PF_condom_msp_1988':   stats.beta_binom(p=.022,n=100),
  'PF_condom_msp_2006':   stats.beta_binom(p=.230,n=483),
  'PF_condom_msp_2016':   stats.beta_binom(p=.416,n= 75),
  'PF_condom_cas_1988':   stats.beta_binom(p=.088,n=316),
  'PF_condom_cas_2006':   stats.beta_binom(p=.598,n=235),
  'PF_condom_cas_2016':   stats.beta_binom(p=.694,n=420),
  'PF_condom_swo_2002':   stats.beta_binom(p=.432,n= 9),
  'PF_condom_swo_2011':   stats.beta_binom(p=.777,n=21),
  'PF_condom_swo_2014':   stats.beta_binom(p=.787,n=14),
  'PF_condom_swr_2002':   stats.beta_binom(p=.337,n=13),
  'PF_condom_swr_2011':   stats.beta_binom(p=.754,n=24),
  'PF_condom_swr_2014':   stats.beta_binom(p=.759,n=11),
  'PF_circum_2050':       stats.beta_binom(p=.724,n=18),
  # beta
  'beta_0':               stats.gamma_p(p=.00075,v=4.2e-8),
  'Rbeta_acute':          stats.gamma_p(p=5.3,v=10),
  'Rbeta_350':            stats.gamma_p(p=1.6,v=.0233),
  'Rbeta_200':            stats.gamma_p(p=8.3,v=4.71),
  'Rbeta_vi_rec':         stats.gamma_p(p=1.45,v=.0658),
  'Rbeta_gud_sus_w':      stats.gamma_p(p=5.3,v=6.8),
  'Rbeta_gud_sus_m:w':    stats.gamma_p(p=1.45,v=.0658),
  'Rbeta_gud_inf':        stats.gamma_p(p=2.9,v=1.45),
  'dur_acute':            stats.gamma_p(p=.1417,v=4.08e-3),
  'P_gud_fsw_l':          stats.beta_binom(p=.16,n=45),
  'RP_gud_fsw_h:l':       stats.gamma_p(p=3,v=.809),
  'RP_gud_2050':          stats.uniform(l=0.2,h=1),
  'iP_gud_cli_fsw:gp':    stats.uniform(l=0,h=1),
  'Rbeta_uvls':           stats.beta_binom(p=.25,n=5),
  # diagnosis
  'dx_2010':              stats.gamma_p(p=.26,v=1.08e-2),
  'aRdx_2020':            stats.gamma_p(p=.30,v=7.09e-2),
  'aRdx_2050':            stats.gamma_p(p=.30,v=7.09e-2),
  'Rdx_mqq':              stats.gamma_p(p=.52,v=2.82e-2),
  'Rdx_cli':              stats.gamma_p(p=.52,v=2.82e-2),
  'Rdx_fsw':              stats.gamma_p(p=2.6,v=1.08),
  # treatment
  'tx_2010':              stats.gamma_p(p=3.0,v=1.67),
  'aRtx_2020':            stats.gamma_p(p=.30,v=7.09e-2),
  'aRtx_2050':            stats.gamma_p(p=.30,v=7.09e-2),
  'Rtx_mqq':              stats.gamma_p(p=.72,v=1.65e-2),
  'Rtx_cli':              stats.gamma_p(p=.72,v=1.65e-2),
  'Rtx_fsw':              stats.gamma_p(p=.72,v=1.65e-2),
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
  NX0   = np.array(243) # REF: WorldBank
  PX0_h = np.array([1,0,0,0,0,0]).reshape([1,1,6,1])
  PX0_c  = np.array([1,0,0,0,0]).reshape([1,1,1,5])
  PX_h_hiv = np.array([0,5,65,30,0,0]).reshape([1,1,6,1])*1e-6 # REF: assume
  PX_h_hiv[:,:,0,:] = 1 - PX_h_hiv.sum()
  P = get_PX_fsw_cli(P)
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
    'X0': NX0 * PX_si[:,:,_,_] * PX0_h * PX0_c,
    'PX_h_hiv': PX_h_hiv,
  }

def get_PX_fsw_cli(P): # [OK]
  P['PX_w']     = .52 # REF: WorldBank
  P['PX_fsw_h'] = .2  # REF: assume
  P['PX_cli_h'] = .2  # REF: assume
  P['F_swo']    = 1   # c.f. get_f()
  P['PX_fsw']   = P['PX_w'] * P['PX_w_fsw']
  P['CF_swo_total'] = P['PX_fsw'] * P['F_swo'] * P['C_swo_fsw_l'] * linear_comb(P['PX_fsw_h'],P['RC_swo_fsw_h:l'],1)
  P['CF_swr_total'] = P['PX_fsw'] * P['F_swr'] * P['C_swr_fsw_l'] * linear_comb(P['PX_fsw_h'],P['RC_swr_fsw_h:l'],1)
  P['PX_cli'] = (P['CF_swo_total'] + P['CF_swr_total']) / P['CF_swq_cli']
  # NOTE: CF_swq_cli is true client average, not CF_swq_cli_l; see get_C()
  return P

def get_birth_death(P): # [OK]
  death = np.array(1/35 + (1-.64)*.0144) # death rate ~ .034 for 15-49 years + non-HIV mortality
  birth = ta.tarray([1980,2000,2010,2020,2050],np.array([.04,.03,.015,.015,P['growth_2050']])+death)
  return {
    'birth_t': birth,
    'death': death,
  }

def get_turnover(P): # [OK]
  t_fsw_l = 1/P['dur_fsw_l'] - P['death']
  t_fsw_h = 1/P['dur_fsw_h'] - P['death']
  t_cli   = 1/P['dur_cli'] - P['death']
  turn = np.zeros((2,4,4)) # s, i.from, i.to
  # fsw -> med & low
  turn[0,2,1] = t_fsw_l * (P['Pturn_fsw_h:m'])
  turn[0,2,0] = t_fsw_l * (1 - P['Pturn_fsw_h:m'])
  turn[0,3,1] = t_fsw_h * (P['Pturn_fsw_h:m'])
  turn[0,3,0] = t_fsw_h * (1 - P['Pturn_fsw_h:m'])
  # clients -> med & low
  Pturn_cli_hm = P['PX_si'][1,1] / P['PX_si'][1,:2].sum()
  turn[1,2,1] = t_cli * (Pturn_cli_hm)
  turn[1,2,0] = t_cli * (1 - Pturn_cli_hm)
  turn[1,3,1] = t_cli * (Pturn_cli_hm)
  turn[1,3,0] = t_cli * (1 - Pturn_cli_hm)
  # balance flows (absolute)
  i = [(3,3,2,2,1),(1,0,1,0,0)]
  turn[:,i[1],i[0]] += turn[:,i[0],i[1]] * P['PX_si'][:,i[0]] / P['PX_si'][:,i[1]]
  ORi = np.exp(P['LORturn_wq_sus:hiv'])
  ORturn_sus_hiv = np.array([[ORi,ORi,1,1],[1,1,1,1]]).reshape((2,4,1))
  return {
    'turn_sii': turn,
    'dur_sii': 1 / (turn.sum(axis=2) + P['death']),
    'ORturn_sus:hiv': ORturn_sus_hiv,
  }

# FOI ----------------------------------------------------------------------------------------------

def get_beta_a(P): # [OK]
  # beta_a dimensions: a,p,s,i,s',i',h',c',(t)
  Rbeta_ar = 10
  Rbeta_as = np.array([[P['Rbeta_vi_rec'],1],[Rbeta_ar,1]]).reshape([2,1,2,1,1,1,1,1])
  Rbeta_as = Rbeta_as / Rbeta_as[0,:].mean()
  P_gud_gp = .07 # REEF: SDHS2006
  P_sexa_l = .6  # REF: (approx) SDHS2006, SHIMS2
  P_gud_cli_h = linear_comb(P['iP_gud_cli_fsw:gp'],  P['P_gud_fsw_l'],P_gud_gp)
  P_gud_cli_l = linear_comb(P['iP_gud_cli_fsw:gp']/2,P['P_gud_fsw_l'],P_gud_gp)
  P_gud = np.array([
    [P_gud_gp*P_sexa_l, P_gud_gp, P['P_gud_fsw_l'], P['P_gud_fsw_l']*P['RP_gud_fsw_h:l']],
    [P_gud_gp*P_sexa_l, P_gud_gp, P_gud_cli_l, P_gud_cli_h],
  ])
  P_gud_t = ta.tarray([1980,2000,2050,2051],[1,1,*2*[P['RP_gud_2050']]])
  Rbeta_gud_sus = np.array([P['Rbeta_gud_sus_w'],P['Rbeta_gud_sus_m:w']]).reshape((2,1))
  Rbeta_h = np.array([0,P['Rbeta_acute'],1,1,P['Rbeta_350'],P['Rbeta_200']]).reshape([1,1,1,1,1,1,6,1])
  Rbeta_c = np.array([1,1,1,P['Rbeta_uvls'],.00]).reshape([1,1,1,1,1,1,1,5])
  return {
    'beta_a': P['beta_0'] * Rbeta_as * Rbeta_h * Rbeta_c,
    'Rbeta_as': Rbeta_as,
    'Rbeta_h': Rbeta_h,
    'Rbeta_gud_sus': Rbeta_gud_sus,
    'P_gud': P_gud,
    'P_gud_t': P_gud_t,
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
  # dimensions: a,p
  F_p   = np.array([ P['F_mcq'], P['F_mcq'], 1, P['F_swr'] ])
  PF_ai = np.array([ P['PF_ai_mcq'],P['PF_ai_mcq'],P['PF_ai_swq'],P['PF_ai_swq'] ])
  F_ap  = F_p.reshape([1,4,1,1,1,1,1,1]) * np.array([1-PF_ai,PF_ai]).reshape([2,4,1,1,1,1,1,1])
  return {
    'F_ap': F_ap,
  }

def check_F(P):
  return (
    P['PF_ai_mcq'] <= P['PF_ai_swq']
  )

def get_C(P): # [OK]
  PX_si = P['PX_si']
  F = np.squeeze(P['F_ap'].sum(axis=0))
  # dimensions: p,s,i
  C_psi = np.zeros((4,2,4))
  # main / spousal
  C_psi[ 0, 0, 0] = .38
  C_psi[ 0, 0, 1] = .32
  C_psi[ 0, 0, 2] = .32
  C_psi[ 0, 0, 3] = .32
  C_psi[ 0, 1, 0] = .35
  C_psi[ 0, 1, 1] = .38
  C_psi[ 0, 1, 2] = .38
  C_psi[ 0, 1, 3] = .38
  # casual
  C_psi[ 1, 0, 0] = .31
  C_psi[ 1, 0, 1] = P['C_cas_am']
  C_psi[ 1, 0, 2] = P['C_cas_am'] * P['RC_cas_fsw:am']
  C_psi[ 1, 0, 3] = P['C_cas_am'] * P['RC_cas_fsw:am']
  C_psi[ 1, 1, 0] = .41
  C_psi[ 1, 1, 1] = P['C_cas_am']
  C_psi[ 1, 1, 2] = P['C_cas_am'] * P['RC_cas_cli:am']
  C_psi[ 1, 1, 3] = P['C_cas_am'] * P['RC_cas_cli:am']
  # FSW
  C_psi[ 2, 0, 2] = P['C_swo_fsw_l']
  C_psi[ 2, 0, 3] = P['C_swo_fsw_l'] * P['RC_swo_fsw_h:l']
  C_psi[ 3, 0, 2] = P['C_swr_fsw_l']
  C_psi[ 3, 0, 3] = P['C_swr_fsw_l'] * P['RC_swr_fsw_h:l']
  # clients
  wPX = (PX_si[1,2] + P['RCF_swq_cli_h:l'] * PX_si[1,3])
  C_psi[ 2, 1, 2] = P['CF_swo_total'] / F[2] / wPX
  C_psi[ 2, 1, 3] = P['CF_swo_total'] / F[2] / wPX * P['RCF_swq_cli_h:l']
  C_psi[ 3, 1, 2] = P['CF_swr_total'] / F[3] / wPX
  C_psi[ 3, 1, 3] = P['CF_swr_total'] / F[3] / wPX * P['RCF_swq_cli_h:l']
  return {
    'C_psi': C_psi,
  }

def get_condom(P): # [OK]
  k = 'PF_condom_' # convenience
  PF_condom_t = ta.tarray([1980,1988,2002,2006,2011,2014,2016,2050],
    [[0,P[k+'msp_1988'],NAN,P[k+'msp_2006'],NAN,NAN,P[k+'msp_2016'],P[k+'msp_2016'] ], # main
     [0,P[k+'cas_1988'],NAN,P[k+'cas_2006'],NAN,NAN,P[k+'cas_2016'],P[k+'cas_2016'] ], # casual
     [0,NAN,P[k+'swo_2002'],NAN,P[k+'swo_2011'],P[k+'swo_2014'],NAN,P[k+'swo_2014'] ], # sw-new
     [0,NAN,P[k+'swr_2002'],NAN,P[k+'swr_2011'],P[k+'swr_2014'],NAN,P[k+'swr_2014'] ]] # sw-reg
  ).reshape([1,4,1,1,1,1,1,1])
  RPF_condom_a = np.array([1,P['RPF_condom_a:v']]).reshape([2,1,1,1,1,1,1,1])
  return {
    'PF_condom_t': PF_condom_t,
    'RPF_condom_a': RPF_condom_a,
    'Rbeta_condom': .26,
  }

def check_condom(P):
  k = 'PF_condom_'
  return (
    # across years (same type)
    P[k+'msp_1988'] < P[k+'msp_2006'] < P[k+'msp_2016'] and
    P[k+'cas_1988'] < P[k+'cas_2006'] < P[k+'cas_2016'] and
    P[k+'swo_2002'] < P[k+'swo_2011'] < P[k+'swo_2014'] and
    P[k+'swr_2002'] < P[k+'swr_2011'] < P[k+'swr_2014'] and
    # across types (same year)
    P[k+'msp_1988'] < P[k+'cas_1988'] and
    P[k+'msp_2006'] < P[k+'cas_2006'] and
    P[k+'msp_2016'] < P[k+'cas_2016'] and
    P[k+'swr_2002'] < P[k+'swo_2002'] and 
    P[k+'swr_2011'] < P[k+'swo_2011'] and 
    P[k+'swr_2014'] < P[k+'swo_2014']
  )

def get_circumcision(P): # [OK]
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
  pref_pii = np.zeros((4,4,4))
  pref_pii[0,0,0] = P['pref_msp_al']
  pref_pii[0,2,2] = P['pref_msp_asw']
  pref_pii[0,2,3] = P['pref_msp_asw']
  pref_pii[0,3,3] = P['pref_msp_asw']
  pref_pii[0,3,2] = P['pref_msp_asw']
  return {
    'pref_pii': pref_pii,
    'mix': np.zeros((4,2,4,2,4)), # initialize
    'mix_mask': np.ones((4,2,4,2,4)),
    't0_tpaf': -1,
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
  # dimensions: s,i,h,c
  # h: sus, acute, >500, <500, <350, <200 (AIDS)
  # c: undiag, diag, unlinked, on art, vls
  dur_h  = np.array([P['dur_acute'],3.5-P['dur_acute'],3.74,5.26]).reshape([1,1,4,1]) # Lodi2011,Mangal2017
  prog_h = 1/dur_h
  unprog_h = np.array([[2.5,.15],[2.5,.15],[2,.15]]).reshape([1,1,3,2]) # Battegay2006,Lawn2006
  death_h = np.array([0,0,.004,.02,.04,.20]).reshape([1,1,6,1]) # Badri2006,Anglaret2012,Mangal2017
  Rdeath_c = np.array([1,1,1,.25,.5]).reshape([1,1,1,5]) # Gabillard2013,Lundgren2015
  return {
    'dur_h': dur_h,
    'prog_h': prog_h,
    'unprog_h': unprog_h,
    'death_hc': death_h * Rdeath_c,
  }

# cascade ------------------------------------------------------------------------------------------

def Rmr(r0,*iRrs,replast=True):
  # helper to simplify sampling monotonic rates from non-zero distributions
  # r0:   base rate
  # iRrs: "incremental" relative rates -- i.e. iRr = 0.2 -> Rr 1.2
  # e.g. Rmr(2, 0.1, 1) -> [2, 2.2, 4.4, 4.4]
  if replast: iRrs = [*iRrs,0]
  Rrs = [1+iRr for iRr in iRrs]
  return np.cumprod([r0,*Rrs])

def get_diag(P):
  # dimensions: s,i,h
  # TODO: reduce sex differences after a time?
  dx_t = ta.tarray([1980,2000,2010,2020,2050,2051],
      [0,0,*Rmr(P['dx_2010'],P['aRdx_2020'],P['aRdx_2050'])])
  Rdx_si = np.array(
    [[           1,           1,P['Rdx_fsw'],P['Rdx_fsw']],
     [P['Rdx_mqq'],P['Rdx_mqq'],P['Rdx_cli'],P['Rdx_cli']]]).reshape([2,4,1])
  Rdx_scen  = np.ones((2,4,1))
  return {
    'dx_t': dx_t,
    'Rdx_si': Rdx_si,
    'Rdx_scen': Rdx_scen,
  }

def get_treat(P):
  # dimensions: s,i,h
  # TODO: reduce sex differences after a time?
  tx_t = ta.tarray([1980,2004,2010,2020,2050,2051],
      [0,0,*Rmr(P['tx_2010'],P['aRtx_2020'],P['aRtx_2050'])])
  Rtx_si = np.array(
    [[           1,           1,P['Rtx_fsw'],P['Rtx_fsw']],
     [P['Rtx_mqq'],P['Rtx_mqq'],P['Rtx_cli'],P['Rtx_cli']]]).reshape([2,4,1])
  Rtx_ht = ta.tarray(
     [1980,2004,2010,2015,2017,2019,2050,2051], [
     [   0,   0,   0,   0,  .1,   1,   1,   1],   # acute:           scale-up 2017-2019
     [   0,   0,   0,   0,  .1,   1,   1,   1],   # cd4 > 500:       scale-up 2017-2019
     [   0,   0,  .1,  .1,   1,   1,   1,   1],   # 350 < cd4 < 500: scale-up 2015-2017
     [   0,   0,   1,   1,   1,   1,   1,   1],   # 200 < cd4 < 350: scale-up 2004-2010
     [   0,   0,   1,   1,   1,   1,   1,   1],   # cd4 < 200:       scale-up 2004-2010
  ]).reshape([1,1,5]) 
  vx   = np.array(1/.36) # median 3 months -> mean ~ 4 months
  unvx_t = ta.tarray([1980,2000,2009,2050,2051],[.2,.2,.1,.04,.04]).reshape([1,1,1]) # NERCHA2014 Figure 15
  retx = np.array(1) # TODO
  Rtx_scen = np.ones((2,4,1))
  Rux_scen = np.ones((2,4,1))
  return {
    'tx_t': tx_t,
    'Rtx_ht': Rtx_ht,
    'Rtx_si': Rtx_si,
    'vx': vx,
    'unvx_t': unvx_t,
    'retx': retx,
    'Rtx_scen': Rtx_scen,
    'Rux_scen': Rux_scen,
  }
