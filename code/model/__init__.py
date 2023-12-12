from model.slicer import slicers

tol = 1e-8

dims = {
  'a': ['V','A'],
  'p': ['M','C','O','R'],
  's': ['W','M'],
  'i': ['L','M','H','Y'],
  'h': ['S','H','1','2','3','A'],
  'c': ['U','D','X','T','V'],
  '*': [''],
}
dimkeys = ''.join(dims.keys())
dimensions = {
  'act':      ['Vaginal','Anal'],
  'partner':  ['Main/Spousal','Casual','Sex Work Occas','Sex Work Regular'],
  'sex':      ['Women','Men'],
  'activity': ['Lowest Activity','Medium Activity','Lower Risk Sex Work','Higher Risk Sex Work'],
  'health':   ['Susceptible','Acute HIV','CD4 > 500','350 < CD4 < 500','200 < CD4 < 350','CD4 < 200 (AIDS)'],
  'care':     ['Undiagnosed','Diagnosed','Virally Unsuppressed','On ART','Virally Suppressed'],
  '*':        [''],
}
