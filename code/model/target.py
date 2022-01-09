import numpy as np
from inspect import signature
from utils import stats
from utils.ops import dictstr,flatten
from model import out

class Target():
  def __init__(self,dist,name,pop,pop2=None,vsop=None,weight=1):
    self.dist = dist
    self.name = name
    self.pop  = pop if not vsop else None
    self.pop1 = pop if vsop else None
    self.pop2 = pop2
    self.vsop = vsop # one of {1/2, 1-2, 1-2/1, 1-2/2}
    self.weight = weight

  def __str__(self):
    return 'Target: {} {{{}}} ~ {}({})'.format(
      self.name,
      self.popstr(),
      self.dist.dist.name,
      dictstr(dict(mu=self.mean(),ci=self.ci()),3),
    )

  def __repr__(self):
    return 'T: {} {{{}}}'.format(
      self.name,
      self.popstr(),
    )

  def popstr(self):
    if not self.vsop:
      return dictstr(self.pop)
    else:
      return out.vs_label(dictstr(self.pop1),dictstr(self.pop2),self.vsop)

  def ll(self,x,interval=None):
    if self.weight:
      if interval:
        xlo,xhi = self.dist.interval(interval)
        return self.weight * (xlo <= x <= xhi)
      else:
        return self.weight * self.dist.logpdf(x)
    else: # False-y
      return self.weight

  def mean(self):
    return self.dist.mean()

  def ci(self,interval=.95):
    return self.dist.interval(interval)

def filter_targets(T,name=None,pop=None,pop1=None,pop2=None,vsop=None):
  eqpop = lambda Tp,p: Tp and all(Tp[k]==v for k,v in p.items() if k in Tp)
  if name is not None:
    T = [Ti for Ti in T if Ti.name in flatten(name)]
  if pop  is not None: T = [Ti for Ti in T if eqpop(Ti.pop, pop )]
  if pop1 is not None: T = [Ti for Ti in T if eqpop(Ti.pop1,pop1)]
  if pop2 is not None: T = [Ti for Ti in T if eqpop(Ti.pop2,pop2)]
  if vsop is not None: T = [Ti for Ti in T if Ti.vsop == vsop]
  return T

def get_model_ll(T,R,t,interval=None):
  ll = 0.0
  for Ti in T:
    if Ti.pop2 is None:
      x = out.by_name(Ti.name)(R,**Ti.pop,tvec=t,aggr=True)
    else:
      x = out.vs_pop(Ti.name,R,pop1=Ti.pop1,pop2=Ti.pop2,vsop=Ti.vsop,tvec=t,aggr=True)
    ll += np.float(Ti.ll(x,interval=interval))
  return ll

def top_q_ll(Rs,top=.1,ll='ll'):
  llcut = np.nanquantile([R[ll] for R in Rs],1-top)
  return [R for R in Rs if R[ll] >= llcut]

def get_all_esw(T=None,**kwds):
  T = T if T is not None else []
  T.extend(get_prevalence_esw())
  T.extend(get_incidence_esw())
  T.extend(get_cd4_esw())
  T.extend(get_cascade_esw())
  T.extend(get_cascade_2020('95-95-95'))
  T.extend(get_pop_total_esw())
  # T.extend(get_prevalence_esw_anc())
  # T.extend(get_spectrum_esw())
  return filter_targets(T,**kwds)

# JK: my analysis
# SR: self-report

def get_prevalence_esw():
  return [
    # Prevalence ratios
    Target(stats.ratio_binom(p1=.61,n1=99,p2=.388,n2=9843),'prevalence',
      dict(t=2011.0,s=0,i=(2,3)),dict(t=2011.0,s=0,i=(0,1,2,3)),vsop='1/2'), # FSW vs W
    Target(stats.gamma_p(p=1.46,v=7.09e-3),'prevalence',
      dict(t=2011.0,s=0,i=3),dict(t=2011.0,s=0,i=2),vsop='1/2'), # FSW.H vs FSW.L
    Target(stats.gamma_p(p=2.3,v=4.52e-2),'prevalence',
      dict(t=2014.0,s=0,i=3),dict(t=2014.0,s=0,i=2),vsop='1/2',weight=0), # FSW.H vs FSW.L
    Target(stats.ratio_binom(p1=.545,n1=373,p2=.382,n2=9412),'prevalence',
      dict(t=2011.0,s=0,i=(1,2,3)),dict(t=2011.0,s=0,i=0),vsop='1/2'), # M.2+ vs M.0-1
    Target(stats.ratio_binom(p1=.281,n1=1515,p2=.232,n2=6733),'prevalence',
      dict(t=2011.0,s=1,i=(1,2,3)),dict(t=2011.0,s=1,i=0),vsop='1/2'), # W.2+ vs W.0-1
    # Baral2014
    Target(stats.beta_binom(p=.61 ,n=   99),'prevalence',dict(t=2011.0,s=0,    i=(2,3))),
    Target(stats.beta_binom(p=.94, n=   65),'prevalence',dict(t=2011.0,s=0,    i=(2))), # JK
    Target(stats.beta_binom(p=.64, n=  260),'prevalence',dict(t=2011.0,s=0,    i=(3))), # JK
    # EswKP2014
    Target(stats.beta_binom(p=.38, n=  781),'prevalence',dict(t=2014.0,s=0,    i=(2,3)),weight=0), # SR
    Target(stats.beta_binom(p=.71, n=  114),'prevalence',dict(t=2014.0,s=0,    i=(2)),  weight=0), # SR,JK
    Target(stats.beta_binom(p=.31, n=  454),'prevalence',dict(t=2014.0,s=0,    i=(3)),  weight=0), # SR,JK
    # EswCOP21
    Target(stats.beta_binom(p=.588,n=   27),'prevalence',dict(t=2021.0,s=0,    i=(2,3))),
    # SDHS2006 (Table B.2) - TODO: include FSW
    Target(stats.beta_binom(p=.311,n= 2688),'prevalence',dict(t=2006.5,s=0,    i=(0,1,2,3))),
    Target(stats.beta_binom(p=.196,n= 1966),'prevalence',dict(t=2006.5,s=1,    i=(0,1,2,3))),
    Target(stats.beta_binom(p=.258,n= 3499),'prevalence',dict(t=2006.5,s=(0,1),i=(0,1,2,3))),
    # SHIMS1 Bicego2013 (Table 3) - TODO: include FSW
    Target(stats.beta_binom(p=.321,n=18172),'prevalence',dict(t=2011.0,s=(0,1),i=(0,1,2,3))),
    Target(stats.beta_binom(p=.388,n= 9843),'prevalence',dict(t=2011.0,s=0,    i=(0,1,2,3))),
    Target(stats.beta_binom(p=.241,n= 8329),'prevalence',dict(t=2011.0,s=1,    i=(0,1,2,3))),
    Target(stats.beta_binom(p=.319,n=16145),'prevalence',dict(t=2011.0,s=(0,1),i=(0))),
    Target(stats.beta_binom(p=.333,n= 1887),'prevalence',dict(t=2011.0,s=(0,1),i=(1,2,3))),
    Target(stats.beta_binom(p=.382,n= 9412),'prevalence',dict(t=2011.0,s=0,    i=(0))),
    Target(stats.beta_binom(p=.545,n=  373),'prevalence',dict(t=2011.0,s=0,    i=(1,2,3))),
    Target(stats.beta_binom(p=.232,n= 6733),'prevalence',dict(t=2011.0,s=1,    i=(0))),
    Target(stats.beta_binom(p=.281,n= 1515),'prevalence',dict(t=2011.0,s=1,    i=(1,2,3))),
    # SHIMS2 (Table 6.3.B) - TODO: include FSW
    Target(stats.beta_binom(p=.272,n= 8533),'prevalence',dict(t=2016.5,s=(0,1),i=(0,1,2,3))),
    Target(stats.beta_binom(p=.343,n= 4878),'prevalence',dict(t=2016.5,s=0,    i=(0,1,2,3))),
    Target(stats.beta_binom(p=.189,n= 3655),'prevalence',dict(t=2016.5,s=1,    i=(0,1,2,3))),
  ]

def get_incidence_esw():
  return [
    # TODO: incidence ratios!!!
    # SHIMS1 Justman2016 (Tables 2,1) - TODO: include FSW
    Target(stats.gamma_p(p=.0310,v=7.95e-6),'incidence',dict(t=2011.0,s=0,    i=(0,1,2,3))),
    Target(stats.gamma_p(p=.0320,v=1.38e-5),'incidence',dict(t=2011.0,s=0,    i=(0))),
    Target(stats.gamma_p(p=.1000,v=1.45e-3),'incidence',dict(t=2011.0,s=0,    i=(1,2,3))),
    Target(stats.gamma_p(p=.0170,v=4.12e-6),'incidence',dict(t=2011.0,s=1,    i=(0,1,2,3))),
    Target(stats.gamma_p(p=.0162,v=1.24e-5),'incidence',dict(t=2011.0,s=1,    i=(0))),
    Target(stats.gamma_p(p=.0380,v=6.45e-5),'incidence',dict(t=2011.0,s=1,    i=(1,2,3))),
    # SHIMS2 (Table 5.3.A) - TODO: include FSW
    Target(stats.gamma_p(p=.0148,v=7.67e-6),'incidence',dict(t=2016.5,s=(0,1),i=(0,1,2,3))),
    Target(stats.gamma_p(p=.0199,v=2.00e-5),'incidence',dict(t=2016.5,s=0,    i=(0,1,2,3))),
    Target(stats.gamma_p(p=.0099,v=1.08e-5),'incidence',dict(t=2016.5,s=1,    i=(0,1,2,3))),
  ]

def get_cd4_esw():
  w = 0
  return [
    # Burtle2012 (Table 1) - omit: not representative sample
    # Target(stats.beta_binom(p=.58,n=200),'Ph',dict(t=2009.2,h=(3,4,5)),weight=w),
    # Target(stats.beta_binom(p=.61,n=771),'Ph',dict(t=2009.4,h=(3,4,5)),weight=w),
    # Target(stats.beta_binom(p=.72,n=200),'Ph',dict(t=2010.2,h=(3,4,5)),weight=w),
    # Target(stats.beta_binom(p=.28,n=200),'Ph',dict(t=2009.2,h=(5)),weight=w),
    # Target(stats.beta_binom(p=.34,n=771),'Ph',dict(t=2009.4,h=(5)),weight=w),
    # Target(stats.beta_binom(p=.45,n=200),'Ph',dict(t=2010.2,h=(5)),weight=w),
    # SHIMS2 (Table 11.3.A)
    Target(stats.beta_binom(p=.440,n=2421),'Ph',dict(t=2016.5,h=(3,4,5)),weight=w),
    Target(stats.beta_binom(p=.077,n=2421),'Ph',dict(t=2016.5,h=(5)),weight=w),
  ]

def get_cascade_esw():
  return [
    # FSW - TODO: review
    Target(stats.beta_binom(p=.776,n= 205),'diagnosed',dict(t=2011.0,s=0,i=(2,3))), # R2P2013
    Target(stats.beta_binom(p=.369,n= 197),'treated_u',dict(t=2011.0,s=0,i=(2,3))), # R2P2013(Adj)
    # Target(stats.beta_binom(p=.596,n= 235),'treated-?',dict(t=2014.0,s=0,i=(2,3))), # EswKP2014(Adj) TODO
    # SHIMS1Tables (Table 6)
    Target(stats.beta_binom(p=.626,n=5742),'diagnosed',dict(t=2011.0,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.691,n=3810),'diagnosed',dict(t=2011.0,s=0,    i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.501,n=1997),'diagnosed',dict(t=2011.0,s=1,    i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.326,n=5742),'treated_u',dict(t=2011.0,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.332,n=3810),'treated_u',dict(t=2011.0,s=0,    i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.314,n=1997),'treated_u',dict(t=2011.0,s=1,    i=(0,1,2,3))), # SR
    # Bicego2013 (SHIMS)
    # Target(stats.beta_binom(p=.692,n=todo),'diagnosed',dict(t=2011.0,s=0,i=(0,1,2,3))), # SR [omit:dup]
    # Target(stats.beta_binom(p=.506,n=todo),'diagnosed',dict(t=2011.0,s=1,i=(0,1,2,3))), # SR [omit:dup]
    # Target(stats.beta_binom(p=.332,n=todo),'treated_u',dict(t=2011.0,s=0,i=(0,1,2,3))), # SR [omit:dup]
    # Target(stats.beta_binom(p=.311,n=todo),'treated_u',dict(t=2011.0,s=1,i=(0,1,2,3))), # SR [omit:dup]
    # EswGARPR2014
    Target(stats.beta_binom(p=.527,n=100),'treated_u',dict(t=2013,s=(0,1),i=(0,1,2,3)),weight=0), # Spectrum
    Target(stats.beta_binom(p=.576,n=100),'treated_u',dict(t=2013,s=0,    i=(0,1,2,3)),weight=0), # Spectrum
    Target(stats.beta_binom(p=.457,n=100),'treated_u',dict(t=2013,s=1,    i=(0,1,2,3)),weight=0), # Spectrum
    # SHIMS2 (Table 10.3.A) - TODO: include FSW
    Target(stats.beta_binom(p=.837,n=2413),'diagnosed',dict(t=2016.5,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.880,n=1687),'diagnosed',dict(t=2016.5,s=0,    i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.744,n= 726),'diagnosed',dict(t=2016.5,s=1,    i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.721,n=2413),'treated_u',dict(t=2016.5,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.757,n=1687),'treated_u',dict(t=2016.5,s=0,    i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.641,n= 726),'treated_u',dict(t=2016.5,s=1,    i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.861,n=2057),'treated_c',dict(t=2016.5,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.861,n=1496),'treated_c',dict(t=2016.5,s=0,    i=(0,1,2,3))), # SR
    Target(stats.beta_binom(p=.862,n= 561),'treated_c',dict(t=2016.5,s=1,    i=(0,1,2,3))), # SR
    # Target(stats.beta_binom(p=.722,n=2413),'treated_u',dict(t=2016.5,s=(0,1),i=(0,1,2,3))), # detect
    # Target(stats.beta_binom(p=.755,n=1686),'treated_u',dict(t=2016.5,s=0,    i=(0,1))), # detect
    # Target(stats.beta_binom(p=.651,n= 727),'treated_u',dict(t=2016.5,s=1,    i=(0,1,2,3))), # detect
    Target(stats.beta_binom(p=.708,n=2423),'vls_u',dict(t=2016.5,s=(0,1),i=(0,1,2,3))),
    Target(stats.beta_binom(p=.748,n=1694),'vls_u',dict(t=2016.5,s=0,    i=(0,1,2,3))),
    Target(stats.beta_binom(p=.623,n= 729),'vls_u',dict(t=2016.5,s=1,    i=(0,1,2,3))),
    Target(stats.beta_binom(p=.909,n=1778),'vls_c',dict(t=2016.5,s=(0,1),i=(0,1,2,3))),
    Target(stats.beta_binom(p=.918,n=1292),'vls_c',dict(t=2016.5,s=0,    i=(0,1,2,3))),
    Target(stats.beta_binom(p=.886,n= 486),'vls_c',dict(t=2016.5,s=1,    i=(0,1,2,3))),
  ]

def get_cascade_2020(which,weight=1):
  def make_targets(p,n,t=2020.0,s=None,i=None):
    if s is None: s = (0,1)
    if i is None: i = (0,1,2,3)
    names = ['diagnosed','treated_c','vls_c']
    pop = dict(t=t,s=s,i=i)
    return [Target(stats.beta_binom(p=p[i],n=n[i]),names[i],pop,weight=weight) for i in range(len(p))]
  if which=='90-90-90':       return make_targets([.900,.900,.900],[4259,4259,4259]) # AIDSinfo
  if which=='95-95-95':       return make_targets([.950,.950,.950],[1650,1650,1650]) # AIDSinfo
  if which=='ssa':            return make_targets([.847,.879,.891],[  31,  32,  32]) # AIDSinfo
  if which=='ssa-lowest':     return make_targets([.168,.465,.733],[ 334, 168,  55]) # AIDSinfo
  if which=='ssa-lowest-fsw': return make_targets([.451,.465,.733],[ 168, 168,  55],s=0,i=(2,3)) # alt-2
# if which=='ssa-lowest-fsw': return make_targets([.168,.465,.733],[ 334, 168,  55],s=0,i=(2,3)) # alt-1
# if which=='ssa-lowest-fsw': return make_targets([.451,.112],     [ 100, 100],     s=0,i=(2,3)) # AIDSinfo
  if which=='90-90-90-fsw':   return make_targets([.900,.900,.900],[4259,4259,4259],s=0,i=(2,3))
  if which=='95-95-95-fsw':   return make_targets([.950,.950,.950],[1650,1650,1650],s=0,i=(2,3))

def get_pop_total_esw():
  return [
    Target(stats.gamma_p(p=Xt/1000,v=Xt/1000),'NX',dict(t=t,s=(0,1),i=(0,1,2,3)),weight=.2)
      for t,Xt in zip(range(1980,2020+1),[
    243151,251090,259122,267588,276979,287513,299011,312027,326013,340165, # 1980-1989
    354047,367323,379398,390945,402993,416073,427619,440611,454193,466912, # 1990-1999
    477968,485089,489755,493025,496291,500406,505239,510583,516321,522181, # 2000-2009
    528097,533686,540179,547254,554459,561694,571423,580164,588576,597650, # 2010-2019
    607854]) if t%5==0 # only every 5th year
  ]

def get_prevalence_esw_anc():
  return [ # UNGASS 2008 (Figure 1) - TODO: include FSW? TODO: find 2010+
    Target(stats.beta_binom(p=Pt,n=1000),'prevalence',dict(t=t,s=0,i=(0,1,2,3)),weight=0)
      for t,Pt in zip(range(1992,2010+1,2),[
    .039,.161,.260,.316,.342,.386,.426,.392,np.nan,.411]) # 1992-2010
  ]

def get_spectrum_esw():
  # from WorldBank
  return [
    Target(stats.beta_binom(p=Pt/100,n=1e4),'prevalence',dict(t=t,s=(0,1),i=(0,1,2,3)),weight=0)
      for t,Pt in zip(range(1990,2020+1),[
      1.1, 2.1, 3.8, 6.3, 9.2,12.3,15.1,17.5,19.3,20.8, # 1990-1999
     21.9,22.9,23.6,24.3,24.8,25.3,25.8,26.2,26.7,27.1, # 2000-2009
     27.6,27.9,28.2,28.6,28.8,28.9,28.8,28.5,28.0,27.4, # 2010-2019
     26.8]) if t%2==0 # only every 2nd year
  ] + [
    Target(stats.beta_binom(p=It/1000,n=1e5),'incidence',dict(t=t,s=(0,1),i=(0,1,2,3)),weight=0)
      for t,It in zip(range(1990,2020+1),[
      6.28,11.60,19.34,28.72,36.85,41.92,43.95,42.80,40.69,37.85, # 1990-1999
     35.40,33.79,32.47,31.41,31.07,30.74,30.35,29.78,29.50,29.28, # 2000-2009
     28.97,27.70,27.46,27.17,26.05,23.30,20.23,16.80,13.48,11.96, # 2010-2019
     10.21]) if t%2==0 # only every 2nd year
  ] + [
    Target(stats.beta_binom(p=Tt/100,n=100),'treated_u',dict(t=t,s=(0,1),i=(0,1,2,3)),weight=0)
      for t,Tt in zip(range(2000,2020+1),[
     0, 1, 2, 2, 4, 7,11,15,20,25,35,39,43,49,59,69,78,87,89,96, # 2000-2020
    98]) if t%2==0 # only every 2nd year
  ] 
  
