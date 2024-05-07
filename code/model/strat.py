# a collection of strata meta-data for plotting

from utils import dict_str

class Strat():
  # object representing a population stratification
  def __init__(self,key,ind,color,label):
    self.key   = key # short unique label
    self.ind   = ind # dimension indices defining the stratification
    self.color = color # color in (r,g,b) format
    self.label = label # long label for plot legends & facets

  def __str__(self):
    return 'Strat: {} {{{}}} [{}]'.format(self.key,dict_str(self.ind),self.label)

  def __repr__(self):
    return 'S: {} {{{}}}'.format(self.key,dict_str(self.ind))

# define and collect all Strats in a dict
strats = {Si.key:Si for Si in [
  Strat('*',     dict(),                   (.40,.00,.40),''),
  Strat('all',   dict(s=(0,1),i=(0,1,2,3)),(.40,.00,.40),'Population Overall'),
  Strat('aq',    dict(s=(0,1),i=(0,1)),    (.60,.20,.60),'Non-Sex Work Women & Men'),
  Strat('w',     dict(s=(0),  i=(0,1,2,3)),(.80,.20,.20),'Women Overall'),
  Strat('wq',    dict(s=(0),  i=(0,1)),    (1.0,.60,.60),'Women (No FSW)'),
  Strat('wl',    dict(s=(0),  i=(0)),      (1.0,.60,.60),'Lowest Activity Women'),
  Strat('wm',    dict(s=(0),  i=(1)),      (.80,.40,.40),'Medium Activity Women'),
  Strat('wh',    dict(s=(0),  i=(1,2,3)),  (.60,.00,.00),'Non-Lowest Activity Women'),
  Strat('fsw.l', dict(s=(0),  i=(2)),      (.60,.20,.20),'Lower Risk FSW'),
  Strat('fsw.h', dict(s=(0),  i=(3)),      (.40,.00,.00),'Higher Risk FSW'),
  Strat('fsw',   dict(s=(0),  i=(2,3)),    (.60,.00,.00),'FSW Overall'),
  Strat('m',     dict(s=(1),  i=(0,1,2,3)),(.20,.20,.80),'Men Overall'),
  Strat('mq',    dict(s=(1),  i=(0,1)),    (.60,.60,1.0),'Men (No Clients)'),
  Strat('ml',    dict(s=(1),  i=(0)),      (.60,.60,1.0),'Lowest Activity Men'),
  Strat('mm',    dict(s=(1),  i=(1)),      (.40,.40,.80),'Medium Activity Men'),
  Strat('mh',    dict(s=(1),  i=(1,2,3)),  (.00,.00,.60),'Non-Lowest Activity Men'),
  Strat('cli.l', dict(s=(1),  i=(2)),      (.20,.20,.60),'Lower Risk Clients'),
  Strat('cli.h', dict(s=(1),  i=(3)),      (.00,.00,.40),'Higher Risk Clients'),
  Strat('cli',   dict(s=(1),  i=(2,3)),    (.00,.00,.60),'Clients Overall'),
  Strat('asw',   dict(s=(0,1),i=(2,3)),    (.60,.00,.60),'FSW & Clients Overall'),
  Strat('asw.l', dict(s=(0,1),i=(2)),      (.60,.20,.60),'Lower Risk FSW & Clients'),
  Strat('asw.h', dict(s=(0,1),i=(3)),      (.40,.00,.40),'Higher Risk FSW & Clients'),
  Strat('ahi',   dict(h=(1)),              (.90,.30,.00),'Acute HIV'),
  Strat('hiv.1', dict(h=(2)),              (.90,.15,.00),'CD4 > 500'),
  Strat('hiv.2', dict(h=(3)),              (.90,.00,.15),'350 < CD4 < 500'),
  Strat('hiv.3', dict(h=(4)),              (.60,.00,.30),'200 < CD4 < 350'),
  Strat('aids',  dict(h=(5)),              (.30,.00,.60),'CD4 < 200 (AIDS)'),
  Strat('>500',  dict(h=(1,2)),            (.90,.15,.00),'CD4 > 500'),
  Strat('<500',  dict(h=(3,4,5)),          (.90,.00,.15),'CD4 < 500'),
  Strat('<350',  dict(h=(4,5)),            (.60,.00,.30),'CD4 < 350'),
  Strat('<200',  dict(h=(5)),              (.30,.00,.60),'CD4 < 200'),
  Strat('udx',   dict(c=(0)),              (.00,.30,.30),'Undiagnosed'),
  Strat('dx',    dict(c=(1)),              (.00,.60,.30),'Diagnosed'),
  Strat('ux',    dict(c=(2)),              (.00,.60,.60),'Unlinked'),
  Strat('tx',    dict(c=(3)),              (.40,.80,.00),'On ART'),
  Strat('vx',    dict(c=(4)),              (.80,.80,.00),'Virally Suppressed'),
  Strat('msp',   dict(p=(0)),              (.26,.04,.41),'Main / Spousal'),
  Strat('cas',   dict(p=(1)),              (.58,.15,.40),'Casual'),
  Strat('swo',   dict(p=(2)),              (.87,.32,.23),'One-Off Sex Work'),
  Strat('swr',   dict(p=(3)),              (.99,.65,.04),'Repeat Sex Work'),
  Strat('swx',   dict(p=(2,3)),            (.93,.48,.13),'Sex Work Overall'),
]}
