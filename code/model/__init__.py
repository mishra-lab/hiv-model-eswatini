from model.slicer import slicers

_ = None

dimkeys = 'apsihc*'
dims = {
  'a': ['V','A'],
  'p': ['L','S','N','R'],
  's': ['W','M'],
  'i': ['L','M','H','Y'],
  'h': ['S','H','1','2','3','A'],
  'c': ['U','D','X','T','V'],
  '*': [''],
}
dimensions = {
  'act':      ['Vaginal','Anal'],
  'partner':  ['Main/Spousal','Casual','New SW','Regular SW'],
  'sex':      ['Women','Men'],
  'activity': ['Lowest Risk','Medium Risk','Lower Risk Sex Work','Higher Risk Sex Work'],
  'health':   ['Susceptible','Acute HIV','CD4 > 500','350 < CD4 < 500','CD4 < 350','AIDS'],
  'care':     ['Undiagnosed','Diagnosed','Lost Contact','Treated','Virally Suppressed'],
  '*':        [''],
}
