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
  Slicer('ALL',   dict(s=(0,1),i=(0,1,2,3)),(.40,.00,.40),'Population Overall'),
  Slicer('AQ',    dict(s=(0,1),i=(0,1)),    (.60,.20,.60),'Non-Sex Work Women & Men'),
  Slicer('W',     dict(s=(0),  i=(0,1,2,3)),(.80,.20,.20),'Women Overall'),
  Slicer('WQ',    dict(s=(0),  i=(0,1)),    (1.0,.60,.60),'Women (No FSW)'),
  Slicer('WL',    dict(s=(0),  i=(0)),      (1.0,.60,.60),'Lowest Risk Women'),
  Slicer('WM',    dict(s=(0),  i=(1)),      (.80,.40,.40),'Medium Risk Women'),
  Slicer('WH',    dict(s=(0),  i=(1,2,3)),  (.60,.00,.00),'Non-Lowest Risk Women'),
  Slicer('FSW.L', dict(s=(0),  i=(2)),      (.60,.20,.20),'Lower Risk FSW'),
  Slicer('FSW.H', dict(s=(0),  i=(3)),      (.40,.00,.00),'Higher Risk FSW'),
  Slicer('FSW',   dict(s=(0),  i=(2,3)),    (.60,.00,.00),'FSW Overall'),
  Slicer('M',     dict(s=(1),  i=(0,1,2,3)),(.20,.20,.80),'Men Overall'),
  Slicer('MQ',    dict(s=(1),  i=(0,1)),    (.60,.60,1.0),'Men (No Clients)'),
  Slicer('ML',    dict(s=(1),  i=(0)),      (.60,.60,1.0),'Lowest Risk Men'),
  Slicer('MM',    dict(s=(1),  i=(1)),      (.40,.40,.80),'Medium Risk Men'),
  Slicer('MH',    dict(s=(1),  i=(1,2,3)),  (.00,.00,.60),'Non-Lowest Risk Men'),
  Slicer('Cli.L', dict(s=(1),  i=(2)),      (.20,.20,.60),'Lower Risk Clients'),
  Slicer('Cli.H', dict(s=(1),  i=(3)),      (.00,.00,.40),'Higher Risk Clients'),
  Slicer('Cli',   dict(s=(1),  i=(2,3)),    (.00,.00,.60),'Clients Overall'),
  Slicer('AHI',   dict(h=(1)),              (.00,.75,.75),'Acute HIV'),
  Slicer('HIV.1', dict(h=(2)),              (.00,.60,.60),'CD4 > 500'),
  Slicer('HIV.2', dict(h=(3)),              (.00,.45,.45),'350 < CD4 < 500'),
  Slicer('HIV.3', dict(h=(4)),              (.00,.30,.30),'200 < CD4 < 350'),
  Slicer('AIDS',  dict(h=(5)),              (.00,.15,.15),'CD4 < 200 (AIDS)'),
  Slicer('>500',  dict(h=(1,2)),            (.00,.60,.60),'CD4 > 500'),
  Slicer('<500',  dict(h=(3,4,5)),          (.00,.45,.45),'CD4 < 500'),
  Slicer('<350',  dict(h=(4,5)),            (.00,.30,.30),'CD4 < 350'),
  Slicer('<200',  dict(h=(5)),              (.00,.15,.15),'CD4 < 200'),
  Slicer('LT',    dict(p=(0)),              (.26,.04,.41),'Long Term'),
  Slicer('ST',    dict(p=(1)),              (.58,.15,.40),'Short Term'),
  Slicer('SWR',   dict(p=(2)),              (.87,.32,.23),'Sex Work Regular'),
  Slicer('SWO',   dict(p=(3)),              (.99,.65,.04),'Sex Work Occasional'),
]}
