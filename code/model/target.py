import numpy as np
from inspect import signature
from utils import deco,stats,flatten,dict_str
from model import out

llmin = -1e6

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
    return 'Target: {} {{{}}} ~ {} @ {}, w = {}'.format(
      self.name,
      self.popstr(),
      self.dist.dist.name,
      dict_str(dict(mu=self.mean(),ci=self.ci())),
      self.weight,
    )

  def __repr__(self):
    return 'T: {} {{{}}}'.format(
      self.name,
      self.popstr(),
    )

  def popstr(self):
    if not self.vsop:
      return dict_str(self.pop)
    else:
      return out.vs_label(dict_str(self.pop1),dict_str(self.pop2),self.vsop)

  @deco.nowarn
  def ll(self,x,interval=None):
    # TODO: describe weights in thesis
    if self.weight:
      if interval:
        xlo,xhi = self.dist.interval(interval)
        return self.weight * (xlo <= x <= xhi)
      else:
        return self.weight * np.maximum(llmin,self.dist.logpdf(x))
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
      x = out.by_name(Ti.name)(R,**Ti.pop,tvec=t)
    else:
      x = out.vs_pop(Ti.name,R,pop1=Ti.pop1,pop2=Ti.pop2,vsop=Ti.vsop,tvec=t)
    ll += float(Ti.ll(x,interval=interval))
  return ll

def top_ll(Rs,top=.1,ll='ll'):
  if isinstance(top,int): top = top / len(Rs)
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
  return filter_targets(T,**kwds)

# (JK): my analysis
# SR: self-report

def get_prevalence_esw():
  return [
    # TODO - check
    # prevalence ratios
    Target(stats.ratio_binom(p1=.605,n1=127,p2=.388,n2=6015),'prevalence',dict(t=2011.0,s=0,i=(2,3)),dict(t=2011.0,s=0,i=(0,1,2,3)),vsop='1/2',weight=10), # FSW vs W
    Target(stats.gamma(m=1.46,sd=.0842),'prevalence',dict(t=2011.0,s=0,i=3),dict(t=2011.0,s=0,i=2),vsop='1/2',weight=10), # FSW.H vs FSW.L JK
    # W.2+ vs W.0-1 & M.2+ vs M.0-1 JK
    Target(stats.invgauss(m=2.01,sd=.129,z=1.76),'prevalence',dict(t=2006.5,s=0,i=(1,2,3)),dict(t=2006.5,s=0,i=0),vsop='1/2'), # SDHS20067 (JK)
    Target(stats.invgauss(m=3.07,sd=.907,z=1.71),'prevalence',dict(t=2006.5,s=1,i=(1,2,3)),dict(t=2006.5,s=1,i=0),vsop='1/2'), # SDHS20067 (JK)
    Target(stats.invgauss(m=1.54,sd=.050,z=1.42),'prevalence',dict(t=2011.0,s=0,i=(1,2,3)),dict(t=2011.0,s=0,i=0),vsop='1/2'), # Bicego2013 (JK)
    Target(stats.invgauss(m=1.25,sd=.038,z=1.18),'prevalence',dict(t=2011.0,s=1,i=(1,2,3)),dict(t=2011.0,s=1,i=0),vsop='1/2'), # Bicego2013 (JK)
    Target(stats.invgauss(m=1.42,sd=.036,z=1.34),'prevalence',dict(t=2016.5,s=0,i=(1,2,3)),dict(t=2016.5,s=0,i=0),vsop='1/2'), # SHIMS2 (JK)
    Target(stats.invgauss(m=1.32,sd=.052,z=1.24),'prevalence',dict(t=2016.5,s=1,i=(1,2,3)),dict(t=2016.5,s=1,i=0),vsop='1/2'), # SHIMS2 (JK)
    # FSW
    Target(stats.betabin(p=.605,n=127),'prevalence',dict(t=2011.0,s=0,i=(2,3)),weight=10), # Baral2014
    Target(stats.betabin(p=.588,n=394),'prevalence',dict(t=2021.0,s=0,i=(2,3)),weight=10), # EswIBBS2022
    # SDHS2006 (Table B.2)
    Target(stats.betabin(p=.258,n=3499),'prevalence',dict(t=2006.5,s=(0,1),i=(0,1,2,3))),
    Target(stats.betabin(p=.311,n=2688),'prevalence',dict(t=2006.5,s=0,    i=(0,1,2,3))),
    Target(stats.betabin(p=.196,n=1976),'prevalence',dict(t=2006.5,s=1,    i=(0,1,2,3))),
    # SHIMS1 Bicego2013 (Table 3) - n adjusted for sampling via '06/16, plus (adj 15-17)
    Target(stats.betabin(p=.280,n=7747),'prevalence',dict(t=2011.0,s=(0,1),i=(0,1,2,3))),
    Target(stats.betabin(p=.342,n=6015),'prevalence',dict(t=2011.0,s=0,    i=(0,1,2,3))),
    Target(stats.betabin(p=.207,n=4977),'prevalence',dict(t=2011.0,s=1,    i=(0,1,2,3))),
    # SHIMS2 (Table C.2)
    Target(stats.betabin(p=.272,n=3620),'prevalence',dict(t=2016.5,s=(0,1),i=(0,1,2,3))),
    Target(stats.betabin(p=.343,n=2994),'prevalence',dict(t=2016.5,s=0,    i=(0,1,2,3))),
    Target(stats.betabin(p=.189,n=2442),'prevalence',dict(t=2016.5,s=1,    i=(0,1,2,3))),
  ]

def get_incidence_esw():
  return [
    # incidence ratioss
    Target(stats.invgauss(m=6.93,sd=3.47,z=4.18),'incidence',dict(t=2011.0,s=0,i=(1,2,3)),dict(t=2011.0,s=0,i=0),vsop='1/2'), # Justman2016
    Target(stats.invgauss(m=4.97,sd=3.39,z=2.85),'incidence',dict(t=2011.0,s=1,i=(1,2,3)),dict(t=2011.0,s=1,i=0),vsop='1/2'), # Justman2016
    # EswIBBS2022 (Table 13)
    Target(stats.invgauss(m=.1195,sd=.0222,z=.0505),'incidence',dict(t=2021.0,s=0,i=(2,3))),
    # SHIMS1 Justman2016 (Tables 2,1) - (adj 15-17)
    Target(stats.skewnorm(m=.0294,sd=.00242,a=  2),'incidence',dict(t=2011.0,s=0,i=(0,1,2,3))),
    Target(stats.skewnorm(m=.0150,sd=.00173,a=  0),'incidence',dict(t=2011.0,s=1,i=(0,1,2,3))),
    Target(stats.skewnorm(m=.0158,sd=.00495,a=-10),'incidence',dict(t=2011.0,s=0,i=0)),
    Target(stats.skewnorm(m=.0075,sd=.00312,a=-10),'incidence',dict(t=2011.0,s=1,i=0)),
    Target(stats.skewnorm(m=.0970,sd=.03619,a=  8),'incidence',dict(t=2011.0,s=0,i=(1,2,3))),
    Target(stats.skewnorm(m=.0342,sd=.00697,a=  2),'incidence',dict(t=2011.0,s=1,i=(1,2,3))),
    # SHIMS2 (Table 5.3.A, 5.3.B)
    Target(stats.skewnorm(m=.0148,sd=.00263),'incidence',dict(t=2016.5,s=(0,1),i=(0,1,2,3)),weight=10),
    Target(stats.skewnorm(m=.0198,sd=.00418),'incidence',dict(t=2016.5,s=0,    i=(0,1,2,3)),weight=10),
    Target(stats.skewnorm(m=.0099,sd=.00306),'incidence',dict(t=2016.5,s=1,    i=(0,1,2,3)),weight=10),
  ]

def get_cd4_esw(w=0):
  return [
    Target(stats.betabin(p=.45, n=3158),'Ph',dict(t=2013.3,h=(4,5)), weight=w), # Jobanputra2015
    # SHIMS2 (Table 11.3.A)
    Target(stats.betabin(p=.440,n=2421),'Ph',dict(t=2016.5,h=(3,4,5)), weight=w),
    Target(stats.betabin(p=.077,n=2421),'Ph',dict(t=2016.5,h=(5)), weight=w),
  ]

def get_cascade_esw():
  return [
    # FSW
    Target(stats.betabin(p=.744,n= 127),'diagnosed',dict(t=2011.0,s=0,i=(2,3))), # R2P2013 JK
    Target(stats.betabin(p=.369,n= 179),'treated_u',dict(t=2011.0,s=0,i=(2,3))), # R2P2013
    # Target(stats.betabin(p=.596,n= 235),'treated-?',dict(t=2014.0,s=0,i=(2,3))), # EswKP2014(Adj) TODO
    Target(stats.betabin(p=.883,n= 411),'diagnosed',dict(t=2021.0,s=0,i=(2,3))), # EswIBBS2022
    Target(stats.betabin(p=.975,n= 363),'treated_c',dict(t=2021.0,s=0,i=(2,3))), # EswIBBS2022
    Target(stats.betabin(p=.861,n= 411),'treated_u',dict(t=2021.0,s=0,i=(2,3))), # EswIBBS2022
    # SHIMS1Tables (Table 6)
    Target(stats.betabin(p=.626,n=5742),'diagnosed',dict(t=2011.0,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.691,n=3810),'diagnosed',dict(t=2011.0,s=0,    i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.501,n=1997),'diagnosed',dict(t=2011.0,s=1,    i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.326,n=5742),'treated_u',dict(t=2011.0,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.332,n=3810),'treated_u',dict(t=2011.0,s=0,    i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.314,n=1997),'treated_u',dict(t=2011.0,s=1,    i=(0,1,2,3))), # SR
    # SHIMS2 (Table 10.3.A)
    Target(stats.betabin(p=.837,n=2413),'diagnosed',dict(t=2016.5,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.880,n=1687),'diagnosed',dict(t=2016.5,s=0,    i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.744,n= 726),'diagnosed',dict(t=2016.5,s=1,    i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.721,n=2413),'treated_u',dict(t=2016.5,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.757,n=1687),'treated_u',dict(t=2016.5,s=0,    i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.641,n= 726),'treated_u',dict(t=2016.5,s=1,    i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.861,n=2057),'treated_c',dict(t=2016.5,s=(0,1),i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.861,n=1496),'treated_c',dict(t=2016.5,s=0,    i=(0,1,2,3))), # SR
    Target(stats.betabin(p=.862,n= 561),'treated_c',dict(t=2016.5,s=1,    i=(0,1,2,3))), # SR
    # Target(stats.betabin(p=.722,n=2413),'treated_u',dict(t=2016.5,s=(0,1),i=(0,1,2,3))), # detect
    # Target(stats.betabin(p=.755,n=1686),'treated_u',dict(t=2016.5,s=0,    i=(0,1))), # detect
    # Target(stats.betabin(p=.651,n= 727),'treated_u',dict(t=2016.5,s=1,    i=(0,1,2,3))), # detect
    Target(stats.betabin(p=.708,n=2423),'vls_u',dict(t=2016.5,s=(0,1),i=(0,1,2,3))),
    Target(stats.betabin(p=.748,n=1694),'vls_u',dict(t=2016.5,s=0,    i=(0,1,2,3))),
    Target(stats.betabin(p=.623,n= 729),'vls_u',dict(t=2016.5,s=1,    i=(0,1,2,3))),
    Target(stats.betabin(p=.909,n=1778),'vls_c',dict(t=2016.5,s=(0,1),i=(0,1,2,3))),
    Target(stats.betabin(p=.918,n=1292),'vls_c',dict(t=2016.5,s=0,    i=(0,1,2,3))),
    Target(stats.betabin(p=.886,n= 486),'vls_c',dict(t=2016.5,s=1,    i=(0,1,2,3))),
  ]

def make_targets_2020(p,n=None,t=2020.0,s=None,i=None,w=1):
  if n is None: n = [100]*len(p)
  if s is None: s = (0,1)
  if i is None: i = (0,1,2,3)
  names = ['diagnosed','treated_c','vls_c']
  pop = dict(t=t,s=s,i=i)
  return [Target(stats.betabin(p=p[i],n=n[i]),names[i],pop,weight=w) for i in range(len(p))]

def get_cascade_2020(which,w=1):
  if which=='90-90-90': return make_targets_2020([.90,.90,.90],w=w)
  if which=='95-95-95': return make_targets_2020([.95,.95,.95],w=w)

def get_pop_total_esw():
  return [ # ages 15-49
    Target(stats.gamma(m=Xt/1000,sd=np.sqrt(Xt/1000)),'NX',dict(t=t,s=(0,1),i=(0,1,2,3)),weight=.1)
      for t,Xt in zip(range(1980,2021+1),[
    243151,251090,259122,267588,276979,287513,299011,312027,326013,340165, # 1980-1989
    354047,367323,379398,390945,402993,416073,427619,440611,454193,466912, # 1990-1999
    477968,485089,489755,493025,496291,500406,505239,510583,516321,522181, # 2000-2009
    528097,533686,540179,547254,554459,561694,571423,580164,588576,597650, # 2010-2019
    607854,618763]) if t%5==0 # only every 5th year
  ]

def get_prevalence_esw_anc(w=0):
  return [
    Target(stats.betabin(p=Pt,n=1000),'prevalence',dict(t=t,s=0,i=(0,1,2,3)),weight=w)
      for t,Pt in zip(# 1992-2010: EswUNGASS2010, 2012-2015: EswHIVPR2015 # TODO: NERCHA2012
    [1992,1994,1996,1998,2000,2002,2004,2006,2008,2010,2012,2013,2014,2015],
    [.039,.161,.260,.316,.342,.386,.426,.392,.381,.411,.344,.384,.350,.370])
  ]

