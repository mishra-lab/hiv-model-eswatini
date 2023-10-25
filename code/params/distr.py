import matplotlib.pyplot as plt
import os
from scipy import optimize
from utils import stats as ss
from utils import fio
import numpy as np

plot = True
debug = False
ttfname = fio.tmpfile('distr-{}.pdf')
p2 = (.025,.975)
p3 = (.025,.5,.975)
p5 = (.025,.25,.5,.75,.975)

# FUNCTIONS ----------------------------------------------------------------------------------------

def pdfmerge():
  os.system(' && '.join([
    'pdftk {} cat output {}'.format(ttfname.format('*'),ttfname.format('all')),
    'mv {} pyplots.pdf'.format(ttfname.format('all')),
    'rm {}'.format(ttfname.format('*')),
  ]))

def section(name):
  print('-'*50+'\n'+name+'\n'+'-'*50)

def wfun(x,w):
  x,w = np.asarray(x),np.asarray(w)
  return np.sum(x*w)/np.sum(w)

def p2r(p,t=1):
  return -np.log(1-p)/t

def q_err(q,distrx,p=p2,pwr=2):
  dq = distrx.ppf(p)
  return np.sum(np.abs(q - dq)**pwr)

def print_distr(distr,x,name,m=None,q=None,p=p2,ok=True):
  xstr = {
    'betabin':     'p = {:.4f}, n = {:5.0f}',
    'gamma':       'm = {:.5f}, sd = {:.5f}',
    'lnorm':       'm = {:.5f}, sd = {:.5f}',
    'skewnorm':    'm = {:.5f}, sd = {:.5f}, a = {:.4f}',
    'invgauss':    'm = {:.5f}, sd = {:.5f}, z = {:.5f}',
    'uniform':     'l = {:.3f}, h = {:.3f}',
    'ratio_binom': 'p1 = {:.4f}, n1 = {:5.0f}, p2 = {:.4f}, n2 = {:5.0f}',
  }[distr.__name__].format(*x)
  qstr = '('+',~'.join(['{:.4f}']*len(p)).format(*distr(*x).ppf(p))+')'
  s = ' --- '.join((distr.__name__,xstr,qstr,name))
  print(s if ok else '(!) FAILED: '+s,flush=True)
  if plot: plot_distr(distr(*x),m,q,p,name)

def fit_distr(distr,m,q,w=1,p=p2,name='',n0=100,sd0=1,a0=0,z0=0,pwr=2):
  m  = wfun(m,w)
  q = [wfun(qi,w) for qi in q]
  dname = distr.__name__
  x0 = {
    'betabin': [m,n0],
    'gamma': [m,sd0],
    'lnorm': [m,sd0],
    'invgauss': [m,sd0,z0],
    'skewnorm': [m,sd0,a0],
  }[distr.__name__]
  xb = {
    'betabin': [(.001,.999),(5,1e5)],
    'gamma': None,
    'lnorm': None,
    'invgauss': [(-1e5,+1e5),(1e-9,1e5),(1e-9,100)],
    'skewnorm': [(-1e5,+1e5),(1e-9,1e5),(-100,100)],
  }[distr.__name__]
  
  jfun = lambda x: q_err(q,distr(*x),p=p,pwr=pwr)
  opt = optimize.minimize(jfun,x0,method='L-BFGS-B',bounds=xb,
    options=dict(ftol=1e-15,gtol=1e-12,eps=1e-9))
  ok = opt.success and np.sqrt(opt.fun) < (1e-2*m)
  print_distr(distr,opt.x,name+('' if ok else ' (!)'),m,q,p,ok)
  if debug and not ok: print(opt)
  return opt.x

def plot_distr(distrx,m=None,q=None,p=p2,name=''):
  dq = distrx.ppf(p); q = dq if q is None else q
  dm = distrx.mean(); m = dm if m is None else m
  xq = (min(*q,*dq),max(*q,*dq))
  x = np.linspace(1.5*xq[0]-.5*xq[1],1.5*xq[1]-.5*xq[0],256)
  d = distrx.pdf(x)
  y = np.quantile(d,.99)*1.01
  plt.figure(figsize=(5,3))
  plt.plot(x,d,color='#CC0033')
  for dk,k,pk in zip((dm,*dq),(m,*q),(False,*p)):
    plt.plot([ k, k],[0,y],'-',color='#0099CC')
    plt.plot([dk,dk],[0,y],':',color='#CC0033')
    if pk: plt.text(k,y*pk,' {:.4f} '.format(distrx.cdf(k)),ha='left',va='center')
  plt.ylim((0,y))
  plt.title(name)
  plt.savefig(ttfname.format(name))
  plt.close()

# PARAMETERS ---------------------------------------------------------------------------------------

def beta():
  section('beta')
  fit_distr(ss.gamma,m=.00072,q=(.0005,.00097),sd0=.0001,name='beta_0') # Boily2009
  fit_distr(ss.gamma,m=.0013,q=(.0005,.0025),sd0=.0005,name='beta_0_adj') # Boily2009
  fit_distr(ss.gamma,m=5.3,q=(1,5.3,15),p=p3,pwr=1,name='Rbeta_acute') # Bellan2015 (CI adj: 15 <- 56)
  fit_distr(ss.gamma,m=.1417,q=(.046,.1417,.5),p=p3,pwr=1,sd0=.1,name='dur_acute') # Bellan2015 - years (CI adj: )
  fit_distr(ss.gamma,m=1.6,q=(1.3,1.9),name='Rbeta_350') # Wawer2005, Boily2009, Donnell2010
  fit_distr(ss.gamma,m=8.3,q=(4.5,13),name='Rbeta_200') # Wawer2005, Boily2009, Donnell2010
  fit_distr(ss.gamma,m=1.45,q=(1,2),name='Rbeta_mtf') # assume
  fit_distr(ss.gamma,m=5.29,q=(1.43,5.29,15),p=p3,name='Rbeta_gud_sus') # Boily2009
  fit_distr(ss.gamma,m=2.9, q=(1.03,2.9,5.69),p=p3,name='Rbeta_gud_inf') # Gray2001
  fit_distr(ss.gamma,m=2.15,q=(0.2,6.0),name='aRbeta_gud_sus_adj') # Boily2009
  fit_distr(ss.gamma,m=0.95,q=(0.2,2.4),name='aRbeta_gud_inf_adj') # Gray2001
  fit_distr(ss.betabin,m=.3,q=(.2,.4),name='P_gud_fsw_l') # (JK) & 18/124 EswIBBS2022
  fit_distr(ss.gamma,m=1.3,q=(1,1.6),name='RP_gud_fsw_h:l') # (JK)
  # @FYI SDHS2006 Table 13.14: ~.07 GUD among wider pop P12M
  # @FYI SDHS2006 Table 13.14: ~.61 GUD among STI/symp/GUD (7% of 11.4%)
  # fit_distr(ss.betabin,m=.25,q=(.015,.65),name='Rbeta_uvls') # Donnell2010 (CI adj)
  return

def circumcision():
  section('circumcision')
  print_distr(ss.betabin,x=[.082,4156],name='circum_2006') # SDHS2006
  print_distr(ss.betabin,x=[.171,8329],name='circum_2011') # SHIMS1 Bicego2013
  print_distr(ss.betabin,x=[.300,3988],name='circum_2016') # SHIMS2
  fit_distr(ss.betabin,m=.7,q=(.5,.9),name='circum_2050') # assume
  return

def condoms():
  section('condoms')
  fit_distr(ss.betabin,m=1-.26,q=(1-.43,1-.13),name='Rbeta_condom') # Giannou2016
  # ai vs vi
  fit_distr(ss.betabin,m=.540,q=(.387,.692),name='PF_condom_ai)') # Owen2020a
  fit_distr(ss.betabin,m=.684,q=(.555,.813),name='PF_condom_vi)') # Owen2020a
  ru = .789; rlse = np.sqrt(1/22 - 1/40 + 1/33 - 1/48); ci = np.exp(np.log(ru)+[-rlse,+rlse])
  fit_distr(ss.betabin,m=.789,q=ci,name='RPF_condom_av') # Owen2020a
  fit_distr(ss.betabin,m=.80,q=(.50,.95),name='RPF_condom_av') # Owen2020a (adj) - TODO
  fit_distr(ss.betabin,m=.80,q=(.7,.9),name='RPF_condom_av (FSW,2011)') # FSW 2011
  fit_distr(ss.betabin,m=.55,q=(.4,.7),name='RPF_condom_av (FSW,2014)') # FSW 2014
  # sex work
  fit_distr(ss.betabin,m=.744*.5,q=(.744*.2,.744),name='PF_condom_swo_2002 (vi)') # SBSS2002
  fit_distr(ss.betabin,m=.600*.5,q=(.600*.2,.600),name='PF_condom_swr_2002 (vi)') # SBSS2002
  fit_distr(ss.betabin,m=.848,q=(.772*.75,.924),name='PF_condom_swo_2011 (vi)') # Baral2014
  fit_distr(ss.betabin,m=.829,q=(.758*.75,.900),name='PF_condom_swr_2011 (vi)') # Baral2014
  # fit_distr(ss.betabin,m=.511,q=(.418,.604), name='PF_condom_npp_2011 (vi)') # Baral2014 [omit]
  print_distr(ss.betabin,x=[.856,100],name='PF_condom_swo_2014') # EswKP2014 (n adj, was 620)
  print_distr(ss.betabin,x=[.885,100],name='PF_condom_swr_2014') # EswKP2014 (n adj, was 595)
  # print_distr(ss.betabin,x=[.806,100],name='PF_condom_npp_2014') # EswKP2014 (n adj, was 395) [omit]
  fit_distr(ss.betabin,m=.884,q=(.549,.95),name='PF_condom_swo_2014 (vi)') # EswKP2014
  fit_distr(ss.betabin,m=.853,q=(.479,.95),name='PF_condom_swr_2014 (vi)') # EswKP2014
  # fit_distr(ss.betabin,m=.801,q=(.647,.90),name='PF_condom_npp_2014 (vi)') # EswKP2014 [omit]
  print_distr(ss.betabin,x=[.500,100],name='PF_condom_swx_2021 (vi)') # EswIBBS2022 always
  print_distr(ss.betabin,x=[.455,11],name='PF_condom_swx_2021 (ai)') # EswIBBS2022 always
  # wider pop
  # fit_distr(ss.betabin,m=.006,q=(.004,.013),name='PF_condom_nsw_wom_1998') # SFHS1988 [omit]
  # fit_distr(ss.betabin,m=.073,q=(.059,.121),name='PF_condom_nsw_men_1998') # SFHS1988 [omit]
  fit_distr(ss.betabin,m=.040,q=(.004,.121),name='PF_condom_cas_1998') # SFHS1988 (JK)
  fit_distr(ss.betabin,m=.600,q=(.535,.660),name='PF_condom_cas_2006') # SDHS2006 Tables 14.7.1, 14.7.2
  fit_distr(ss.betabin,m=.693,q=(.649,.737),name='PF_condom_cas_2016') # SHIMS2 Tables Table 15.4.A, Table 15.4.B
  fit_distr(ss.betabin,m=.222,q=(.194,.269),name='PF_condom_msp_2006') # SDHS2006 Tables 14.7.1, 14.7.2
  fit_distr(ss.betabin,m=.414,q=(.308,.529),name='PF_condom_msp_2016') # SHIMS2 Tables Table 15.4.A, Table 15.4.B
  return

def PX():
  section('pop sizes')
  # FSW
  fit_distr(ss.betabin,m=.029, q=(.007,.065),name='PX_w_fsw_2011') # EswKP2014 Table 4
  fit_distr(ss.betabin,m=.0243,q=(.0117,.0502),name='PX_w_fsw_2021') # EswKP2014 Table 4
  fit_distr(ss.betabin,m=.164,q=(.096,.278),name='PX_w_h') # (JK)
  # fit_distr(ss.betabin,m=.252,q=(.150,.440),name='PX_m_h') # (JK) [omit]
  fit_distr(ss.betabin,m=.130,q=(.100,.170),name='PX_m_m') # (JK)
  # fit_distr(ss.betabin,m=.130,q=(.021,.385),name='PX_m_cli') # (JK) [omit]
  return

def CF():
  section('partnership numbers')
  # fit_distr(ss.betabin,m=.38,q=(.21,.57),name='C_msp_wl') # (JK) [omit]
  # fit_distr(ss.betabin,m=.32,q=(.09,.55),name='C_msp_wh') # (JK) [omit]
  # fit_distr(ss.betabin,m=.35,q=(.23,.50),name='C_msp_ml') # (JK) [omit]
  fit_distr(ss.betabin,m=.37,q=(.25,.50),name='C_msp_xl') # (JK) synthesis
  fit_distr(ss.betabin,m=.35,q=(.20,.55),name='C_cas_xl') # (JK) synthesis
  fit_distr(ss.gamma,m=1.5,q=(1.2,2),name='C_cas_wm') # assume
  print_distr(ss.uniform,x=[.25,1],name='RC_cas_cli:wm') # assume
  fit_distr(ss.gamma,m=8.4,q=(6,11),name='CF_swr_l') # (JK)
  print_distr(ss.gamma,x=[3.5,.70],name='C_swo_fsw_l') # (JK) Baral2014,EswKP2014
  print_distr(ss.gamma,x=[6.0,1.2],name='C_swr_fsw_l') # (JK) Baral2014,EswKP2014
  print_distr(ss.gamma,x=[14,2.8], name='C_swo_fsw_h') # (JK) Baral2014,EswKP2014
  print_distr(ss.gamma,x=[21,4.2], name='C_swr_fsw_h') # (JK) Baral2014,EswKP2014
  # fit_distr(ss.gamma,m=4.1,q=(2.5,6.0),name='C_swo_fsw_l') # (JK) Baral2014,EswKP2014 [omit]
  # fit_distr(ss.gamma,m=2.0,q=(1.6,2.5),name='RC_swo_fsw_h:l') # (JK) Baral2014,EswKP2014 [omit]
  # fit_distr(ss.gamma,m=8.4,q=(6.0,11.0),name='C_swr_fsw_l') # (JK) Baral2014,EswKP2014 [omit]
  # fit_distr(ss.gamma,m=8.4*5/6,q=(6.0*5/6,11.0*5/6),name='C_swr_fsw_l_adj') # (JK) Baral2014,EswKP2014 [omit]
  # fit_distr(ss.gamma,m=1.5,q=(1.3,1.7),name='RC_swr_fsw_h:l') # (JK) Baral2014,EswKP2014 [omit]
  fit_distr(ss.gamma,m=36,q=(18,72),name='KF_swx_cli') # assume
  fit_distr(ss.gamma,m=2.0,q=(1.6,2.5),name='RKF_swx_cli_h:l') # assume
  section('sex frequency')
  fit_distr(ss.gamma,m=52*1.5,q=(.5*52,3*52),name='F_msp') # Shisana2005,Delva2013
  print_distr(ss.uniform,x=[12,48],name='F_swr') # assume
  fit_distr(ss.gamma,m=.05,q=(.006,.165),sd0=.01,name='PF_ai_mcx') # Owen2017
  print_distr(ss.uniform,x=[1,2],name='RPF_ai_swx:mcx') # Owen2017
  # fit_distr(ss.gamma,m=.08,q=(.006,.292),sd0=.01,name='PF_ai_swx_old') # Owen2017 [omit]
  # fit_distr(ss.gamma,m=.08,q=(.024,.159),sd0=.01,name='PF_ai_swx') # Owen2020 [omit]
  return

def mixing():
  section('mixing')
  # fit_distr(ss.gamma,m=9,q=(2,19),name='pref_mcx_swx') # assume [omit]
  # fit_distr(ss.gamma,m=2,q=(1.5,3),name='pref_msp_xl') # (JK) [omit]
  print_distr(ss.uniform,x=[-2,+2],name='lpref_mxc_swx') # assume
  print_distr(ss.uniform,x=[-2,+2],name='lpref_msp_xl') # assume
  return

def dur():
  section('durations: risk groups')
  # fit_distr(ss.betabin,m=.383,q=(.275,.383,.491),p=p3,name='dur_fsw_0-2')  # Baral2014 [prelim]
  # fit_distr(ss.betabin,m=.321,q=(.236,.321,.407),p=p3,name='dur_fsw_3-5')  # Baral2014 [prelim]
  # fit_distr(ss.betabin,m=.202,q=(.132,.202,.271),p=p3,name='dur_fsw_6-10') # Baral2014 [prelim]
  # fit_distr(ss.betabin,m=.094,q=(.044,.094,.144),p=p3,name='dur_fsw_11+')  # Baral2014 [prelim]
  fit_distr(ss.gamma,m=4.06,q=(2.29,6.34),name='dur_fsw') # (JK) - Knight2024* bias
  fit_distr(ss.gamma,m=0.5,q=(1/6,1.00),name='dur_sw_h') # assume
  # fit_distr(ss.invgauss,m=2,q=(1.54,2,3.25),p=p3,name='Rdur_fsw_hl') # (JK)
  # fit_distr(ss.gamma,m=4.07,q=(2.96, 5.48),name='dur_fsw_lr') # (JK)
  # fit_distr(ss.gamma,m=9.33,q=(6.30,13.13),name='dur_fsw_hr') # (JK)
  # fit_distr(ss.betabin,m=.5,q=(.4,.6),name='Pturn_e_fsw') # assume
  # fit_distr(ss.betabin,.075,q=(.056,.075,.100),p=p3,name='ybs_ever_15-24') # Hodgins2022 [prelim]
  # fit_distr(ss.betabin,.088,q=(.065,.088,.117),p=p3,name='ybs_ever_25-34') # Hodgins2022 [prelim]
  # fit_distr(ss.betabin,.077,q=(.057,.077,.103),p=p3,name='ybs_ever_35-54') # Hodgins2022 [prelim]
  # fit_distr(ss.betabin,.051,q=(.036,.051,.071),p=p3,name='ybs_p12m_15-24') # Hodgins2022 [prelim]
  # fit_distr(ss.betabin,.039,q=(.027,.039,.056),p=p3,name='ybs_p12m_25-34') # Hodgins2022 [prelim]
  # fit_distr(ss.betabin,.022,q=(.015,.022,.032),p=p3,name='ybs_p12m_35-54') # Hodgins2022 [prelim]
  fit_distr(ss.gamma,m=8,q=(4,15),name='dur_cli') # (JK)
  fit_distr(ss.gamma,m=0.5,q=(1/6,1.00),name='dur_cli_h') # assume
  fit_distr(ss.gamma,m=.2,q=(.05,.50),sd0=.1,name='turn_xm_xl') # assume
  fit_distr(ss.betabin,m=.75,q=(.50,.90),name='Pturn_fsw_m:l') # assume
  fit_distr(ss.betabin,m=.60,q=(.25,.90),name='Pturn_cli_m:l') # assume
  section('durations: partnership types')
  print_distr(ss.uniform,x=[14.5,18.5],name='dur_msp') # (JK)
  # print_distr(ss.uniform,x=[.01,.02],name='pcr_msp') # (JK) [omit]
  fit_distr(ss.gamma,m=.75,q=(.25,1.5),name='dur_cas') # (JK)
  fit_distr(ss.gamma,m=.5,q=(1/6,1),name='dur_swr') # (JK)
  return

def dx():
  section('HIV test p12m')
  # FSW
  fit_distr(ss.gamma,m=p2r(.617),q=(p2r(.556),p2r(.675)),name='dx_fsw_2011') # Baral2014
  print_distr(ss.gamma,[p2r(.617)/p2r(.468)-1,.2],name='aRdx_fsw:w_2011_adj') # (JK)
  print_distr(ss.betabin,[.746,100],name='test_fsw_2014') # EswKP2014 (adj n was 781)
  fit_distr(ss.gamma,m=p2r(.75),q=(p2r(.657),p2r(.826)),name='dx_fsw_2014_adj') # EswKP2014 (adj)
  print_distr(ss.gamma,[p2r(.746)/p2r(.571)-1,.2],name='aRdx_fsw:w_2016_adj') # (JK)
  # 2000 / 2002
  print_distr(ss.betabin,[.133,5642],name='test_ever_women_2000') # EswMICS2000
  fit_distr(ss.betabin,m=.1,q=(.05,.15),name='dx_w_2002') # assume
  # SDHS2006
  fit_distr(ss.betabin,m=p2r(.219),q=(p2r(.206),p2r(.233)),name='dx_w_2006') # SDHS2006
  fit_distr(ss.betabin,m=p2r(.089),q=(p2r(.078),p2r(.100)),name='dx_m_2006') # SDHS2006
  print_distr(ss.betabin,[p2r(.219),100],name='dx_w_2006_adj') # SDHS2006 (adj n)
  print_distr(ss.gamma,[p2r(.089)/p2r(.219),.1],name='Rdx_m:w_2006_adj') # SDHS2006 (adj n)
  # SHIMS
  print_distr(ss.betabin,[.5051,100],name='test_w_2011') # SHIMS1T (adj n was 9838)
  print_distr(ss.betabin,[.3167,100],name='test_m_2011') # SHIMS1T (adj n was 8318)
  print_distr(ss.betabin,[.468,100],name='test_w_2011_adj') # SHIMS1T (adj 15-17 & n was 9838)
  print_distr(ss.betabin,[.284,100],name='test_m_2011_adj') # SHIMS1T (adj 15-17 & n was 8318)
  fit_distr(ss.gamma,m=p2r(.468),q=(p2r(.372),p2r(.566)),sd0=.01,name='dx_w_2011_adj') # (JK)
  fit_distr(ss.gamma,m=p2r(.284),q=(p2r(.200),p2r(.376)),sd0=.01,name='dx_m_2011_adj') # (JK)
  print_distr(ss.gamma,[p2r(.284)/p2r(.468),.1],name='Rdx_m:w_2011_adj') # (JK)
  # 2016
  print_distr(ss.betabin,[.529,100],name='test_all_2016') # SHIMS2 (adj n was 9075) [omit]
  print_distr(ss.betabin,[p2r(.571),100],name='dx_w_2016') # SHIMS2 (adj n was 5127)
  print_distr(ss.betabin,[p2r(.478),100],name='dx_m_2016') # SHIMS2 (adj n was 3948)
  print_distr(ss.gamma,[p2r(.571)/p2r(.505)-1,.05],name='Rdx_w_16:11') # (JK)
  print_distr(ss.gamma,[p2r(.478)/p2r(.317)-1,.05],name='Rdx_m_16:11') # (JK) [omit]
  print_distr(ss.gamma,[p2r(.479)/p2r(.571),.1],name='Rdx_m:w_2016_adj') # (JK)
  return

def tx():
  section('ART initiation rate amond diagnosed')
  fit_distr(ss.gamma,m=1.5,q=(.5,3),name='tx_2010') # assume
  fit_distr(ss.gamma,m=9,q=(6,12),name='tx_2012') # assume
  fit_distr(ss.gamma,m=.75,q=(.5,1),name='revx_2010') # assume
  fit_distr(ss.gamma,m=.6,q=(.33,1),sd0=.1,name='ivx') # assume
  return

# TARGETS ------------------------------------------------------------------------------------------

def prevalence():
  section('prevalence')
  # prevalence ratios
  print_distr(ss.ratio_binom,x=[.605,127,.280,7747],p=p3,name='prev_fsw_w_2011') # Baral2014, Bicego2013
  fit_distr(ss.gamma,m=1.46,q=(1.30,1.63),name='prev_fsw_hl_2011') # Baral2014 (JK)
  fit_distr(ss.gamma,m=2.3, q=(1.92,2.75),name='prev_fsw_hl_2014') # EswKP2014 (JK) [omit] SR
  args = dict(p=p5,pwr=2,sd0=.1)
  fit_distr(ss.invgauss,m=1.979,q=(1.843,1.918,1.977,2.063,2.335),**args,name='prev_women_hl_2006_adj') # SDHS2006 (adj wp)
  fit_distr(ss.invgauss,m=2.571,q=(2.155,2.361,2.572,2.927,5.280),**args,name='prev_men_hl_2006_adj')   # SDHS2006 (adj wp)
  fit_distr(ss.invgauss,m=1.526,q=(1.467,1.501,1.526,1.558,1.660),**args,name='prev_women_hl_2011_adj') # Bicego2013 (adj wp & 15-17)
  fit_distr(ss.invgauss,m=1.236,q=(1.200,1.220,1.236,1.261,1.344),**args,name='prev_men_hl_2011_adj')   # Bicego2013 (adj wp & 15-17)
  fit_distr(ss.invgauss,m=1.410,q=(1.368,1.391,1.409,1.434,1.506),**args,name='prev_women_hl_2016_adj') # SHIMS2 (adj wp)
  fit_distr(ss.invgauss,m=1.305,q=(1.258,1.283,1.305,1.337,1.454),**args,name='prev_men_hl_2016_adj')   # SHIMS2 (adj wp)
  # FSW
  fit_distr(ss.betabin,m=.605,q=(.521,.690),name='prev_fsw_2011') # Baral2014 (RDS-adj)
  fit_distr(ss.betabin,m=.588,q=(.539,.636),name='prev_fsw_2021') # EswIBBS2022 (RDS-adj)
  # lr: 0-1 pp6m; nlr: 2+ p6m
  # SDHS2006 (Tables B.2, 14.3, 14.6, 14.7)
  fit_distr(ss.betabin,m=.259,q=(.244,.273),name='prev_all_2006')
  fit_distr(ss.betabin,m=.311,q=(.294,.329),name='prev_women_2006')
  fit_distr(ss.betabin,m=.197,q=(.179,.214),name='prev_men_2006')
  print_distr(ss.betabin,x=[wfun([.033,.298,.330],[1971,922,4696]),1971+922+4696],name='prev_all_lr_2006') # (unadj wp)
  print_distr(ss.betabin,x=[wfun([.052,.388,.358],[ 780,577,2989]), 780+577+2989],name='prev_women_lr_2006') # (unadj wp
  print_distr(ss.betabin,x=[wfun([.021,.149,.280],[1191,345,1707]),1191+345+1707],name='prev_men_lr_2006') # (unadj wp)
  fit_distr(ss.betabin,m=0.268,q=(.227,.287),name='prev_women_lr_2006_adj') # (adj wp) [omit->pr]
  fit_distr(ss.betabin,m=0.141,q=(.065,.167),name='prev_men_lr_2006_adj')   # (adj wp) [omit->pr]
  print_distr(ss.betabin,x=[wfun([.380,.414],[524,55]),524+55],name='prev_all_nlr_2006')   # [omit->pr]
  print_distr(ss.betabin,x=[wfun([.523,.645],[ 68, 4]), 68+ 4],name='prev_women_nlr_2006') # [omit->pr]
  print_distr(ss.betabin,x=[wfun([.359,.380],[456,50]),456+50],name='prev_men_nlr_2006')   # [omit->pr]
  # SHIMS1 Bicego2013 (Table 3) - adj n by fratio '06 & '16
  print_distr(ss.betabin,x=[.321,18172],name='prev_all_2011')       # (unadj n & 15-17)
  print_distr(ss.betabin,x=[.388, 9843],name='prev_women_2011')     # (unadj n & 15-17)
  print_distr(ss.betabin,x=[.241, 8329],name='prev_men_2011')       # (unadj n & 15-17)
  print_distr(ss.betabin,x=[.280, 7747],name='prev_all_2011_adj')   # (adj n & 15-17)
  print_distr(ss.betabin,x=[.342, 6015],name='prev_women_2011_adj') # (adj n & 15-17)
  print_distr(ss.betabin,x=[.207, 4977],name='prev_men_2011_adj')   # (adj n & 15-17)
  print_distr(ss.betabin,x=[wfun([.198,.360],[4062,12083]),4062+12083],name='prev_all_lr_2011') # (unadj wp)
  print_distr(ss.betabin,x=[wfun([.338,.392],[1794,7618]),1794+7618],name='prev_women_lr_2011') # (unadj wp)
  print_distr(ss.betabin,x=[wfun([.087,.306],[2267,4466]),2267+4466],name='prev_men_lr_2011')   # (unadj wp)
  fit_distr(ss.betabin,m=.315,q=(.289,.328),name='prev_women_lr_2011_adj') # (adj wp & 15-17) [omit->pr]
  fit_distr(ss.betabin,m=.195,q=(.180,.201),name='prev_men_lr_2011_adj')   # (adj wp & 15-17) [omit->pr]
  print_distr(ss.betabin,x=[.333,1887],name='prev_all_nlr_2011')       # (unadj n & 15-17)
  print_distr(ss.betabin,x=[.545, 373],name='prev_women_nlr_2011')     # (unadj n & 15-17)
  print_distr(ss.betabin,x=[.281,1515],name='prev_men_nlr_2011')       # (unadj n & 15-17)
  print_distr(ss.betabin,x=[.290, 792],name='prev_all_nlr_2011_adj')   # (adj n & 15-17) [omit->pr]
  print_distr(ss.betabin,x=[.481, 216],name='prev_women_nlr_2011_adj') # (adj n & 15-17) [omit->pr]
  print_distr(ss.betabin,x=[.241, 955],name='prev_men_nlr_2011_adj')   # (adj n & 15-17) [omit->pr]
  # SHIMS 2 (Table C.2, [6.3.B], 15.3.A)
  fit_distr(ss.betabin,m=.272,q=(.258,.287),name='prev_all_2016')
  fit_distr(ss.betabin,m=.343,q=(.326,.360),name='prev_women_2016')
  fit_distr(ss.betabin,m=.189,q=(.173,.204),name='prev_men_2016')
  print_distr(ss.betabin,x=[wfun([.324,.321],[1723,6164]),1723+6164],name='prev_all_lr_2016') # (unadj wp)
  print_distr(ss.betabin,x=[wfun([.361,.367],[1355,3848]),1355+3848],name='prev_women_lr_2016') # (unadj wp)
  print_distr(ss.betabin,x=[wfun([.226,.255],[ 368,2316]), 368+2316],name='prev_men_lr_2016') # (unadj wp)
  fit_distr(ss.betabin,m=.321,q=(.300,.331),name='prev_women_lr_2016_adj') # (adj wp) [omit->pr]
  fit_distr(ss.betabin,m=.175,q=(.157,.181),name='prev_men_lr_2016_adj')   # (adj wp) [omit->pr]
  print_distr(ss.betabin,x=[.287,914],name='prev_all_nlr_2016')   # [omit->pr]
  print_distr(ss.betabin,x=[.453,263],name='prev_women_nlr_2016') # [omit->pr]
  print_distr(ss.betabin,x=[.228,651],name='prev_men_nlr_2016')   # [omit->pr]
  return

def incidence():
  section('incidence')
  # EswIBBS2022 (Table 13)
  q = ss.betabin(p=-np.log(1-30/676)/141*365,n=676).ppf(p3) * [141/160,1,141/119]
  fit_distr(ss.invgauss,m=q[1],q=q,p=p3,sd0=.1,name='inc_fsw_2021')
  # incidence ratios
  args = dict(p=(.05,.25,.50,.75,.95)) # p != p5
  fit_distr(ss.invgauss,m=5.891,q=(4.398,5.088,5.891,7.310,13.394),**args,name='inc_women_hl_2011_adj') # Bicego2013 (adj wp & 15-17)
  fit_distr(ss.invgauss,m=3.939,q=(2.908,3.395,3.939,5.055,10.864),**args,name='inc_men_hl_2011_adj')   # Bicego2013 (adj wp & 15-17)
  # SHIMS1 Justman2016 (Tables 2,1)
  fit_distr(ss.skewnorm,m=.031,q=(.026,.037),a0=2,name='inc_women_2011')        # (unadj 15-17)
  fit_distr(ss.skewnorm,m=.017,q=(.013,.021),name='inc_men_2011')               # (unadj 15-17)
  fit_distr(ss.skewnorm,m=.0294,q=(.0252,.0347),a0=2,name='inc_women_2011_adj') # (adj 15-17)
  fit_distr(ss.skewnorm,m=.0150,q=(.0116,.0184),name='inc_men_2011_adj')        # (adj 15-17)
  fit_distr(ss.skewnorm,m=[.011,.036],q=([.006,.030],[.022,.044]),w=[789,4135],a0= 3,name='inc_women_lr_2011') # (unadj wp & 15-17)
  fit_distr(ss.skewnorm,m=[.004,.020],q=([.001,.015],[.018,.027]),w=[909,2946],a0=10,name='inc_men_lr_2011')  # (unadj wp & 15-17)
  fit_distr(ss.skewnorm,m=.0163,q=(.0040,.0224),a0=-10,name='inc_women_lr_2011_adj')  # (adj wp & 15-17) [omit->pr]
  fit_distr(ss.skewnorm,m=.0084,q=(.0001,.0117),a0=-10,name='inc_men_lr_2011_adj')    # (adj wp & 15-17) [omit->pr]
  fit_distr(ss.skewnorm,m=.100,q=(.050,.192),a0=8,name='inc_women_nlr_2011')          # (unadj 15-17)
  fit_distr(ss.skewnorm,m=.038,q=(.025,.056),a0=2,name='inc_men_nlr_2011')            # (unadj 15-17)
  fit_distr(ss.skewnorm,m=.0948,q=(.0476,.1829),a0=8,name='inc_women_nlr_2011_adj')   # (adj 15-17) [omit->ir]
  fit_distr(ss.skewnorm,m=.0339,q=(.0221,.0494),a0=2,name='inc_men_nlr_2011_adj')     # (adj 15-17) [omit->ir]
  # SHIMS2 (Table C.1)
  fit_distr(ss.skewnorm,m=.0148,q=(.0096,.0199),name='inc_all_2016')
  fit_distr(ss.skewnorm,m=.0199,q=(.0116,.0280),name='inc_women_2016')
  fit_distr(ss.skewnorm,m=.0099,q=(.0039,.0159),name='inc_men_2016')
  return

def diagnosed():
  section('diagnosed')
  # EswR2P2013 (Table 10)
  fit_distr(ss.betabin,m=.450,q=(.395,.506),name='diag_u_fsw_2011') # n.b. among FSW not FSW-LHIV + prev (simulate) -> below
  fit_distr(ss.skewnorm,m=.741,q=(.617,.741,.898),p=p3,pwr=1,a0=2,name='diag_fsw_2011')
  # EswIBBS2022 (Table 13)
  print_distr(ss.betabin,x=[363/411,411],name='diag_fsw_2021')
  # SHIMS1T (Table 6)
  print_distr(ss.betabin,x=[.626,5807],name='diag_all_2011')
  print_distr(ss.betabin,x=[.691,3810],name='diag_women_2011')
  print_distr(ss.betabin,x=[.501,1997],name='diag_men_2011')
  # SHIMS2 (Table C.7)
  fit_distr(ss.betabin,m=.861,q=(.847,.876),name='diag_all_2016')
  fit_distr(ss.betabin,m=.902,q=(.886,.918),name='diag_women_2016')
  fit_distr(ss.betabin,m=.773,q=(.740,.806),name='diag_men_2016')
  return

def treated():
  section('treated')
  # EswKP2014 (Table 10)
  fit_distr(ss.betabin,m=.369,q=(.301,.442),name='treat_c_fsw_2011')
  fit_distr(ss.skewnorm,m=.274,q=(.209,.274,.357),p=p3,pwr=1,a0=2,name='treat_u_fsw_2011')
  print_distr(ss.betabin,x=[354/363,363],name='treat_c_fsw_2021')
  print_distr(ss.betabin,x=[354/411,411],name='treat_u_fsw_2021')
  # SHIMS1T (Table 6)
  print_distr(ss.betabin,x=[.521,3635],name='treat_c_all_2011')
  print_distr(ss.betabin,x=[.480,2633],name='treat_c_women_2011')
  print_distr(ss.betabin,x=[.627,1002],name='treat_c_men_2011')
  print_distr(ss.betabin,x=[.319,5807],name='treat_u_all_2011')
  print_distr(ss.betabin,x=[.332,3810],name='treat_u_women_2011')
  print_distr(ss.betabin,x=[.314,1997],name='treat_u_men_2011')
  # SHIMS2 (Tables C.6 C.7)
  fit_distr(ss.betabin,m=.878,q=(.860,.896),name='treat_c_all_2016')
  fit_distr(ss.betabin,m=.875,q=(.854,.896),name='treat_c_women_2016')
  fit_distr(ss.betabin,m=.884,q=(.852,.916),name='treat_c_men_2016')
  fit_distr(ss.betabin,m=.756,q=(.736,.775),name='treat_u_all_2016')
  fit_distr(ss.betabin,m=.789,q=(.768,.811),name='treat_u_women_2016')
  fit_distr(ss.betabin,m=.683,q=(.647,.720),name='treat_u_men_2016')
  return

def vls():
  section('vls')
  # SHIMS2 (Tables C.6 C.7)
  fit_distr(ss.betabin,m=.903,q=(.890,.916),name='vls_c_all_2016')
  fit_distr(ss.betabin,m=.914,q=(.899,.928),name='vls_c_women_2016')
  fit_distr(ss.betabin,m=.876,q=(.844,.909),name='vls_c_men_2016')
  fit_distr(ss.betabin,m=.682,q=(.661,.704),name='vls_u_all_2016')
  fit_distr(ss.betabin,m=.721,q=(.697,.745),name='vls_u_women_2016')
  fit_distr(ss.betabin,m=.599,q=(.561,.637),name='vls_u_men_2016')
  return

def cascade_ideal():
  section('cascade: 90-90-90 / 95-95-95')
  # assume
  fit_distr(ss.betabin,m=.90**1,q=(.99*.90**1,1.01*.90**1),name='90^1')
  fit_distr(ss.betabin,m=.90**2,q=(.99*.90**2,1.01*.90**2),name='90^2')
  fit_distr(ss.betabin,m=.90**3,q=(.99*.90**3,1.01*.90**3),name='90^3')
  fit_distr(ss.betabin,m=.95**1,q=(.99*.95**1,1.01*.95**1),name='95^1')
  fit_distr(ss.betabin,m=.95**2,q=(.99*.95**2,1.01*.95**2),name='95^2')
  fit_distr(ss.betabin,m=.95**3,q=(.99*.95**3,1.01*.95**3),name='95^3')
  return

# MAIN ---------------------------------------------------------------------------------------------

# beta()
# circumcision()
# condoms()
# PX()
# CF()
# mixing()
# dur()
# dx()
# tx()
# prevalence()
# incidence()
# diagnosed()
# treated()
# vls()
# cascade_ideal()

if plot: pdfmerge()