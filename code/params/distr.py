import matplotlib.pyplot as plt
from scipy import optimize,stats
from utils import stats as ss
import numpy as np

# FUNCTIONS ----------------------------------------------------------------------------------------

def section(name):
  print('-'*50+'\n'+name+'\n'+'-'*50)

def wfun(x,w):
  x,w = np.asarray(x),np.asarray(w)
  return np.sum(x*w)/np.sum(w)

def ci_fun(dist,interval=.95):
  k = (1-interval)/2
  return dist.ppf([k,1-k])

def ci_err(ci,dist,interval=.95):
  k = (1-interval)/2
  return np.sqrt(np.sum((ci-ci_fun(dist,interval))**2))

def dist_str(dist,x,name,success=True):
  if success:
    if dist == 'beta_binom':
      print('[beta_binom]: p = {:.3f}, n = {:5.0f} --- {}'.format(x[0],x[1],name),flush=True)
    if dist == 'gamma':
      print('[gamma]: p = {:.4f}, v = {:.2e} --- {}'.format(x[0],x[1],name),flush=True)
    if dist == 'lnorm':
      print('[lnorm]: m = {:.4f}, s = {:.2e} --- {}'.format(x[0],x[1],name),flush=True)
  else:
    print('{}: (!) failed to converge --- {} @ {:.4f}, {:.4f}'.format(name,dist,*x),flush=True)

def plot_dist(dist,ci=None,bound=(0,1)):
  if ci is None: ci = ci_fun(dist)
  mm = (2*ci[0]-ci[1],2*ci[1]-ci[0])
  u = np.linspace(max(bound[0],mm[0]),min(bound[1],mm[1]),1024) if bound else\
      np.linspace(mm[0],mm[1],1024)
  # u = np.arange(max(0,1-2*(1-ci[0])),min(1,2*ci[1]),1/1024)
  plt.plot(u,dist.pdf(u),'b')
  plt.plot([ci[0],ci[0]],[0,dist.pdf(ci[0])],'r')
  plt.plot([ci[1],ci[1]],[0,dist.pdf(ci[1])],'r')
  p = [dist.cdf(ci[0]),1-dist.cdf(ci[1])]
  plt.text(ci[0],2*dist.pdf(ci[0]),'{:.4f}'.format(p[0]),ha='right')
  plt.text(ci[1],2*dist.pdf(ci[1]),'{:.4f}'.format(p[1]),ha='left')
  plt.text(ci[0]/2+ci[1]/2,0,'{:.4f}'.format(1-p[0]-p[1]),ha='center')
  plt.show()

def find_beta_binom(m,ci,w=1,n0=100,name='',interval=.95,plot=False):
  m    = wfun(m,w) # unused in optimize
  ci   = [wfun(ci[0],w),wfun(ci[1],w)]
  x0   = [m,n0]
  jfun = lambda x: ci_err(ci,ss.beta_binom(p=x[0],n=x[1]),interval=interval)
  opt  = optimize.minimize(jfun,x0,method='L-BFGS-B',bounds=[(.002,.99),(5,1e5)])
  opt.success = opt.success and opt.fun < 1e-3
  dist_str('beta_binom',opt.x,name,opt.success)
  if plot:
    plot_dist(ss.beta_binom(p=opt.x[0],n=opt.x[1]),ci)

def find_gamma(m,ci,w=1,log10v=-5,name='',interval=.95,plot=False):
  m    = wfun(m,w)
  ci   = [wfun(ci[0],w),wfun(ci[1],w)]
  x0   = log10v
  jfun = lambda x: ci_err(ci,ss.gamma_p(p=m,v=10**x),interval=interval)
  opt  = optimize.minimize(jfun,x0,method='L-BFGS-B')
  opt.success = opt.success and opt.fun < m/50
  dist_str('gamma',[m,10**opt.x[0]],name,opt.success)
  if plot:
    plot_dist(ss.gamma_p(p=m,v=10**opt.x[0]),ci,bound=False)

# SECTIONS -----------------------------------------------------------------------------------------

def beta(plot=False):
  section('beta')
  find_gamma(m=.00072,ci=(.0005,.001),name='beta_0',plot=plot) # Boily2009
  find_gamma(m=5.3,ci=(1,13),name='Rbeta_acute',plot=plot) # Bellan2015 (CI adj)
  find_gamma(m=.1417,ci=(.045,.3),name='dur_acute',plot=plot) # Bellan2015 (CI adj) - years
  # find_gamma(m=1.7,ci=(.55,3.5),name='Rbeta_acute',plot=plot) # Bellan2015 (CI adj) - months
  find_gamma(m=1.6,ci=(1.3,1.9),name='Rbeta_350',plot=plot) # Wawer2005, Boily2009, Donnell2010
  find_gamma(m=8.3,ci=(4.5,13),name='Rbeta_200',plot=plot) # Wawer2005, Boily2009, Donnell2010
  find_gamma(m=1.45,ci=(1,2),name='Rbeta_mtf',plot=plot) # assume
  find_gamma(m=7.6,ci=(1.4,19.5), name='Rbeta_gud_s',plot=plot) # Boily2009
  find_gamma(m=2.9,ci=(1.03,5.69),name='Rbeta_gud_i',plot=plot) # Gray2001
  find_beta_binom(m=.2,ci=(.1,.4),name='P_gud_fsw_l',plot=plot) # params/fsw
  find_gamma(m=3,ci=(1.5,5),name='P_gud_fsw_hl',plot=plot) # params/fsw
  # @FYI SDHS2006 Table 13.14: ~.07 GUD among wider pop P12M
  # @FYI SDHS2006 Table 13.14: ~.61 GUD among STI/symp/GUD (7% of 11.4%)
  find_beta_binom(m=.25,ci=(.015,.65),name='Rbeta_uvls',plot=plot) # Donnell2010 (CI adj)

def A(plot=False):
  section('sex acts')
  find_gamma(6,(1,15),name='CA_swq',plot=True)
  exit()
  find_beta_binom(.1,(.006,.165),name='PA_ai_mc: Wider Pop ai/vi',plot=plot) # Owen2017
  find_beta_binom(.1,(.006,.292),name='PA_ai_sw: Higher-Risk ai/vi',plot=plot) # Owen2017
  find_gamma(78,(26,156),name='A_mc: sex freq',plot=plot) # Delva2013
  find_beta_binom(6/12,(1/12,9/12),name='dur_cas (?)',plot=plot) # REF: TODO
  # find_gamma(.56,(.25,1),name='(?)',plot=plot) # ???

def condoms(plot=True):
  section('condoms')
  # ai vs vi
  find_beta_binom(.540,(.387,.692),name='PA_condom_ai (FSW)',plot=plot) # Owen2020a
  find_beta_binom(.684,(.555,.813),name='PA_condom_vi (FSW)',plot=plot) # Owen2020a
  ru = .789; rlse = np.sqrt(1/22 - 1/40 + 1/33 - 1/48); ci = np.exp(np.log(ru)+[-rlse,+rlse])
  find_beta_binom(.789,(.664,.938),name='RPA_condom_av (FSW)',plot=plot) # Owen2020a
  find_beta_binom(.80,(.50,.95),name='RPA_condom_av',plot=plot) # Owen2020a (adj)
  find_beta_binom(.80,(.7,.9),name='RPA_condom_av (FSW 2011)',plot=plot) # FSW 2011
  find_beta_binom(.55,(.4,.7),name='RPA_condom_av (FSW 2014)',plot=plot) # FSW 2014
  # sex work
  find_beta_binom(.744*.5,(.744*.2,.744),name='PA_cond_new (2002,vi)') # SBSS2002
  find_beta_binom(.600*.5,(.600*.2,.600),name='PA_cond_reg (2002,vi)') # SBSS2002
  find_beta_binom(.848,(.772*.75,.924),name='PA_cond_new (2011,vi)') # Baral2014
  find_beta_binom(.835,(.758*.75,.900),name='PA_cond_reg (2011,vi)') # Baral2014 (p adj, was .829)
  find_beta_binom(.51,(.418,.604), name='PA_cond_npp (2011,vi)') # Baral2014 (p adj, was .511) [omit]
  print(ci_fun(ss.beta_binom(p=.856,n=100))) # EswKP2014 PA_condom_new' (n adj, was 620)
  print(ci_fun(ss.beta_binom(p=.885,n=100))) # EswKP2014 PA_condom_reg' (n adj, was 595)
  print(ci_fun(ss.beta_binom(p=.806,n=100))) # EswKP2014 PA_condom_npp' (n adj, was 395) [omit]
  find_beta_binom(.8,(.549,.95),name='PA_cond_new (2014,vi)',plot=plot) # EswKP2014 (p adj, was .884)
  find_beta_binom(.8,(.479,.95),name='PA_cond_reg (2014,vi)',plot=plot) # EswKP2014 (p adj, was .853)
  find_beta_binom(.8,(.647,.90),name='PA_cond_npp (2014,vi)',plot=plot) # EswKP2014 (p adj, was .801) [omit]
  # wider pop
  find_beta_binom(.023,(.004,.059),name='PA_cond_msp (1998)',plot=plot) # SFHS1988 (p adj, was .006)
  find_beta_binom(.088,(.059,.121),name='PA_cond_cas (1998)',plot=plot) # SFHS1988 (p adj, was .073)
  find_beta_binom(.60, (.535,.660),name='PA_cond_cas (2006)',plot=plot) # SDHS2006 Tables 14.7.1, 14.7.2
  find_beta_binom(.693,(.649,.737),name='PA_cond_cas (2016)',plot=plot) # SHIMS2 Tables Table 15.4.A, Table 15.4.B
  find_beta_binom(.222,(.194,.269),name='PA_cond_msp (2006)',plot=plot) # SDHS2006 Tables 14.7.1, 14.7.2
  find_beta_binom(.414,(.308,.529),name='PA_cond_msp (2016)',plot=plot) # SHIMS2 Tables Table 15.4.A, Table 15.4.B

def cascade_theo():
  section('cascade: 90-90-90 / 95-95-95')
  # assume
  find_beta_binom(m=.900,ci=(.891,.909),name='90^1')
  find_beta_binom(m=.810,ci=(.802,.818),name='90^2')
  find_beta_binom(m=.729,ci=(.722,.736),name='90^3')
  find_beta_binom(m=.950,ci=(.939,.960),name='95^1')
  find_beta_binom(m=.902,ci=(.893,.911),name='95^2')
  find_beta_binom(m=.857,ci=(.849,.866),name='95^3')

def cascade_ssa():
  section('cascade: SSA')
  # AIDSinfo
  find_beta_binom(m=.89,ci=(.72,.95),name='ESA 2020 Diagnosed')
  find_beta_binom(m=.77,ci=(.60,.92),name='ESA 2020 Treated (u)')
  find_beta_binom(m=.70,ci=(.57,.83),name='ESA 2020 VLS (u)')
  find_beta_binom(m=.77,ci=(.63,.94),name='WCA 2020 Diagnosed')
  find_beta_binom(m=.73,ci=(.58,.90),name='WCA 2020 Treated (u)')
  find_beta_binom(m=.59,ci=(.49,.72),name='WCA 2020 VLS (u)')
  find_beta_binom(m=[.89,.77],ci=([.72,.63],[.95,.94]),w=[20.6,4.7],                  name='SSA 2020 Diagnosed')
  find_beta_binom(m=[.77,.73],ci=([.60,.58],[.92,.90]),w=[20.6,4.7],                  name='SSA 2020 Treated (u)')
  dist_str('beta_binom',[wfun([.77/.89,.73/.77],[18.3,3.6]),wfun([34,23],[18.3,3.6])],name='SSA 2020 Treated (c)')
  find_beta_binom(m=[.70,.59],ci=([.57,.49],[.83,.72]),w=[20.6,4.7],                  name='SSA 2020 VLS (u)')
  dist_str('beta_binom',[wfun([.70/.77,.59/.73],[16.0,3.5]),wfun([34,23],[16.0,3.5])],name='SSA 2020 VLS (c)')

def incidence_esw():
  section('incidence')
  # SHIMS1 Justman2016 Tables 2,1
  find_gamma(m=.031,       ci=(.026,.037),                           name='2011 Women Overall')
  find_gamma(m=[.011,.036],ci=([.006,.030],[.022,.044]),w=[789,4135],name='2011 Women 0-1 PP6M')
  find_gamma(m=.100,       ci=(.050,.192),                           name='2011 Women 2+  PP6M')
  find_gamma(m=.017,       ci=(.013,.021),                           name='2011 Men Overall')
  find_gamma(m=[.004,.020],ci=([.001,.015],[.018,.027]),w=[909,2946],name='2011 Men 0-1 PP6M')
  find_gamma(m=.038,       ci=(.025,.056),                           name='2011 Men 2+  PP6M')
  # SHIMS2 Table 5.3.A
  find_gamma(m=.0148,ci=(.0093,.0203),name='2016 Overall')
  find_gamma(m=.0199,ci=(.0109,.0288),name='2016 Women Overall')
  find_gamma(m=.0099,ci=(.0032,.0166),name='2016 Men Overall')

def prevalence_esw(plot=False):
  section('prevalence')
  # prevalence ratios
  find_gamma(m=1.46,ci=(1.30,1.63),name='2011 FSW.H / FSW.L',plot=plot)
  find_gamma(m=2.3,ci=(1.92,2.75),name='2011 FSW.H / FSW.L',plot=plot)
  exit()
  # FSW
  find_beta_binom(m=.61,ci=(.514,.705),name='2011 FSW',plot=plot) # Baral2014
  # find_beta_binom(.605,(.521,.690),name='FSW 2011') # ??? [same]
  find_beta_binom(m=.588,ci=(.40,.76),name='2021 FSW',plot=plot) # EswCOP21 Table 2.1.1; CI assume
  find_beta_binom(m=.311,ci=(.294,.329), name='2006 Women Overall',plot=plot) # SDHS2006
  find_beta_binom(m=.197,ci=(.179,.2141),name='2006 Men Overall',  plot=plot) # SDHS2006
  find_beta_binom(m=.259,ci=(.244,.273), name='2006 All Overall',  plot=plot) # SDHS2006
  # SHIMS1 Bicego2013 (Table 3)
  dist_str('beta_binom',[.388, 9843],                               '2011 Overall')
  dist_str('beta_binom',[.241, 8329],                               '2011 Women Overall')
  dist_str('beta_binom',[.321,18172],                               '2011 Men Overall')
  dist_str('beta_binom',[wfun([.198,.360],[4062,12083]),4062+12083],'2011 All 0-1 PP6M')
  dist_str('beta_binom',[.333, 1887],                               '2011 All 2+  PP6M')
  dist_str('beta_binom',[wfun([.338,.392],[1794,7618]),1794+7618],  '2011 Women 0-1 PP6M')
  dist_str('beta_binom',[.545,  373],                               '2011 Women 2+  PP6M')
  dist_str('beta_binom',[wfun([.087,.306],[2267,4466]),2267+4466],  '2011 Men 0-1 PP6M')
  dist_str('beta_binom',[.281, 1515],                               '2011 Men 2+  PP6M')
  # SHIMS 2 (Table 6.3.B)
  dist_str('beta_binom',[.272,8533],'2016 Overall')
  dist_str('beta_binom',[.343,4878],'2016 Women Overall')
  dist_str('beta_binom',[.189,3655],'2016 Men Overall')

def circumcision():
  section('circumcision')
  dist_str('beta_binom',[.171,8329],'2011') # SHIMS1 Bicego2013
  dist_str('beta_binom',[.300,3988],'2016') # SHIMS2
  find_beta_binom(.6,(.5,.9),name='circum 2050',plot=plot) # assume

def cascade_rates(plot=True):
  find_gamma(m=.26,ci=(.10,.50), name='dx_2010 (wqq)',plot=plot) # assume
  find_gamma(m=.52,ci=(.25,.90), name='Rdx_mqq',plot=plot) # assume
  find_gamma(m=.72,ci=(.50,1.0), name='Rdx_cli',plot=plot) # assume
  find_gamma(m=2.6,ci=(1.0,5.0), name='Rdx_fsw',plot=plot) # assume
  find_gamma(m=3.0,ci=(1.0,6.0), name='tx',plot=plot) # assume
  find_gamma(m=.16,ci=(.01,.50), name='1+Rtx 2017/2010',plot=plot) # assume
  find_gamma(m=.45,ci=(.25,.70), name='tx',plot=plot) # assume [OLD]
  find_gamma(m=.30,ci=(.01,1.0), name='1+Rdx 2017/2010',plot=plot) # assume [OLD]

# MAIN ---------------------------------------------------------------------------------------------


find_beta_binom(.16,(.13,.21),name='dx 2020 ssa-lowest',plot=True)
find_beta_binom(.46,(.39,.54),name='tx 2020 ssa-lowest',plot=True)
find_beta_binom(.72,(.61,.84),name='tx 2020 ssa-lowest',plot=True)

# A()
# condoms()
# exit()
# find_gamma(m= 4.9,ci=(4.0, 5.8),name='dur_fsw_l',plot=True)
# find_gamma(m=10.0,ci=(9.0,11.0),name='dur_fsw_h',plot=True)
# find_gamma(m=10,ci=(6,15.0),name='dur_cli',plot=True)
# find_beta_binom(m=.7,ci=(.5,.9),name='Pturn_h:m',plot=True)
# prevalence_esw()
# cascade_rates()

# TODO
# def PX(plot=False):
  # find_beta_binom(m=.7,ci=(.5,.95),name='turnover FSW -> med',plot=True)
  # find_beta_binom(m=.029,ci=(.006,.065),name='PX_fsw',n0=300,plot=True) # EswKP2014 Table 4
  # find_beta_binom(m=.07,ci=(.02,.15),name='PX_*m (JK)',plot=True)
  # find_beta_binom(m=.1,ci=(.023,.229),name='PX_*m (SDHS2006)',plot=True)
  # find_beta_binom(m=.1,ci=(.038,.182),name='PX_*m (SHIMS1)',plot=True) # Table 3
  # find_beta_binom(m=.1,ci=(.048,.195),name='PX_*m (SHIMS2)',plot=True) # Table 15.3.A
  # find_beta_binom(m=.,ci=(.,.),         name='C_ms (SDHS2006)',plot=True) # Table 14.7.1-2
  # find_beta_binom(m=.665,ci=(.492,.765),name='C_ms (SHIMS1)',plot=True) # Table 3
  # find_beta_binom(m=.720,ci=(.648,.841),name='C_ms (SHIMS2)',plot=True) # Table 15.3.A
  # find_gamma(m=4.1,ci=(2.5,6),name='C_new_fswl',  plot=True)
  # find_gamma(m=8.4,ci=(6,11), name='C_reg_fswl',  plot=True)
  # find_gamma(m=2.0,ci=(1.6,2.5),name='C_new_fsw_hl',plot=True)
  # find_gamma(m=1.5,ci=(1.3,1.7), name='C_reg_fsw_hl',plot=True)
  # find_gamma(m=2,ci=(1.2,3), name='C_reg_cli',plot=True)
  # find_gamma(m=2,ci=(1.2,3), name='A_reg',plot=True)
  # find_gamma(m=3.5,ci=(2,6),name='DX_fsw')
  # find_gamma(m=2.8,ci=(1.2,5), name='C_reg_cli',plot=True)
# find_gamma(m=2.8,ci=(1.2,5), name='A_swq_cli',plot=True)
