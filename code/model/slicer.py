from utils import dict_str

class Slicer():

  def __init__(self,key,pop,color,label):
    self.key   = key
    self.pop   = pop
    self.color = color
    self.label = label

  def __str__(self):
    return 'Slicer: {} {{{}}} [{}]'.format(self.key,dict_str(self.pop),self.label)

  def __repr__(self):
    return 'S: {} {{{}}}'.format(self.key,dict_str(self.pop))

slicers = {Si.key:Si for Si in [
  Slicer('*',     dict(),                   (.40,.00,.40),''),
  Slicer('all',   dict(s=(0,1),i=(0,1,2,3)),(.40,.00,.40),'Population Overall'),
  Slicer('aq',    dict(s=(0,1),i=(0,1)),    (.60,.20,.60),'Non-Sex Work Women & Men'),
  Slicer('w',     dict(s=(0),  i=(0,1,2,3)),(.80,.20,.20),'Women Overall'),
  Slicer('wq',    dict(s=(0),  i=(0,1)),    (1.0,.60,.60),'Women (No FSW)'),
  Slicer('wl',    dict(s=(0),  i=(0)),      (1.0,.60,.60),'Lowest Activity Women'),
  Slicer('wm',    dict(s=(0),  i=(1)),      (.80,.40,.40),'Medium Activity Women'),
  Slicer('wh',    dict(s=(0),  i=(1,2,3)),  (.60,.00,.00),'Non-Lowest Activity Women'),
  Slicer('fsw.l', dict(s=(0),  i=(2)),      (.60,.20,.20),'Lower Risk FSW'),
  Slicer('fsw.h', dict(s=(0),  i=(3)),      (.40,.00,.00),'Higher Risk FSW'),
  Slicer('fsw',   dict(s=(0),  i=(2,3)),    (.60,.00,.00),'FSW Overall'),
  Slicer('m',     dict(s=(1),  i=(0,1,2,3)),(.20,.20,.80),'Men Overall'),
  Slicer('mq',    dict(s=(1),  i=(0,1)),    (.60,.60,1.0),'Men (No Clients)'),
  Slicer('ml',    dict(s=(1),  i=(0)),      (.60,.60,1.0),'Lowest Activity Men'),
  Slicer('mm',    dict(s=(1),  i=(1)),      (.40,.40,.80),'Medium Activity Men'),
  Slicer('mh',    dict(s=(1),  i=(1,2,3)),  (.00,.00,.60),'Non-Lowest Activity Men'),
  Slicer('cli.l', dict(s=(1),  i=(2)),      (.20,.20,.60),'Lower Risk Clients'),
  Slicer('cli.h', dict(s=(1),  i=(3)),      (.00,.00,.40),'Higher Risk Clients'),
  Slicer('cli',   dict(s=(1),  i=(2,3)),    (.00,.00,.60),'Clients Overall'),
  Slicer('asw',   dict(s=(0,1),i=(2,3)),    (.60,.00,.60),'FSW & Clients Overall'),
  Slicer('asw.l', dict(s=(0,1),i=(2)),      (.60,.20,.60),'Lower Risk FSW & Clients'),
  Slicer('asw.h', dict(s=(0,1),i=(3)),      (.40,.00,.40),'Higher Risk FSW & Clients'),
  Slicer('ahi',   dict(h=(1)),              (.90,.30,.00),'Acute HIV'),
  Slicer('hiv.1', dict(h=(2)),              (.90,.15,.00),'CD4 > 500'),
  Slicer('hiv.2', dict(h=(3)),              (.90,.00,.15),'350 < CD4 < 500'),
  Slicer('hiv.3', dict(h=(4)),              (.60,.00,.30),'200 < CD4 < 350'),
  Slicer('aids',  dict(h=(5)),              (.30,.00,.60),'CD4 < 200 (AIDS)'),
  Slicer('>500',  dict(h=(1,2)),            (.90,.15,.00),'CD4 > 500'),
  Slicer('<500',  dict(h=(3,4,5)),          (.90,.00,.15),'CD4 < 500'),
  Slicer('<350',  dict(h=(4,5)),            (.60,.00,.30),'CD4 < 350'),
  Slicer('<200',  dict(h=(5)),              (.30,.00,.60),'CD4 < 200'),
  Slicer('udx',   dict(c=(0)),              (.00,.30,.30),'Undiagnosed'),
  Slicer('dx',    dict(c=(1)),              (.00,.60,.30),'Diagnosed'),
  Slicer('ux',    dict(c=(2)),              (.00,.60,.60),'Unlinked'),
  Slicer('tx',    dict(c=(3)),              (.40,.80,.00),'On ART'),
  Slicer('vx',    dict(c=(4)),              (.80,.80,.00),'Virally Suppressed'),
  Slicer('msp',   dict(p=(0)),              (.26,.04,.41),'Main / Spousal'),
  Slicer('cas',   dict(p=(1)),              (.58,.15,.40),'Casual'),
  Slicer('swo',   dict(p=(2)),              (.87,.32,.23),'Sex Work Occasional'),
  Slicer('swr',   dict(p=(3)),              (.99,.65,.04),'Sex Work Regular'),
]}
