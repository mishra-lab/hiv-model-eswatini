source('utils/ops.r')
source('utils/plot.r')

uid = '2022-02-20'
N = list(cal=100000,batchs=10,sens=10,topfit=.01)
pops = list(
  'wl'    = list(clr=rgb(1.,.6,.6),lab='Women Low'),
  'wm'    = list(clr=rgb(.8,.4,.4),lab='Women Med'),
  'fsw.l' = list(clr=rgb(.6,.2,.2),lab='LR FSW'),
  'fsw.h' = list(clr=rgb(.4,.0,.0),lab='HR FSW'),
  'ml'    = list(clr=rgb(.6,.6,1.),lab='Men Low'),
  'mm'    = list(clr=rgb(.4,.4,.8),lab='Men Med'),
  'cli.l' = list(clr=rgb(.2,.2,.6),lab='LR Clients'),
  'cli.h' = list(clr=rgb(.0,.0,.4),lab='HR Clients'))
parts = list(
  'msp' = list(clr=rgb(.26,.04,.41),lab='Main / Spousal'),
  'cas' = list(clr=rgb(.58,.15,.40),lab='Casual'),
  'swo' = list(clr=rgb(.87,.32,.23),lab='Occasional Sex Work'),
  'swr' = list(clr=rgb(.99,.65,.04),lab='Regular Sex Work'))
groups = list(
  'all' = list(clr=rgb(.4,.4,.4),lab='Overall'),
  'aq'  = list(clr=rgb(.9,.7,.0),lab='Lower Risk'),
  'fsw' = list(clr=rgb(.6,.0,.0),lab='FSW'),
  'cli' = list(clr=rgb(.0,.0,.6),lab='Clients'))
steps = list(
  'diagnosed' = list(clr=rgb(.0,.0,.0),lab='Diagnosed among HIV+'),
  'treated.c' = list(clr=rgb(.0,.0,.0),lab='Treated among diagnosed'),
  'vls.c'     = list(clr=rgb(.0,.0,.0),lab='VLS among treated'),
  'treated.u' = list(clr=rgb(.0,.0,.0),lab='Treated among HIV+'),
  'vls.u'     = list(clr=rgb(.0,.0,.0),lab='VLS among HIV+'))
cases = list(
  'fsw-cli+' = list(clr=rgb(.9,.3,.3),id='-+',lab='Left Behind: FSW'),
  'fsw+cli-' = list(clr=rgb(.3,.3,.9),id='+-',lab='Left Behind: Clients'),
  'fsw-cli-' = list(clr=rgb(.8,.3,.8),id='--',lab='Left Behind: FSW & Clients'),
  'fsw+cli+' = list(clr=rgb(.9,.7,.0),id='++',lab='Left Behind: Neither'),
  'base'     = list(clr=rgb(.4,.4,.4),id='bc',lab='Base Case'))
cases.sens = list(
  'sens' = list(id='rl',lab='Random Lower'),
  'base' = list(id='bc',lab='Base Case')
)
cget = function(defs,key){
  return(unname(sapply(defs,function(def){def[[key]]})))
}
case.ids = cget(cases,'id')
clr = list(
  from   = cget(pops,'clr'),
  to     = cget(pops,'clr'),
  part   = cget(parts,'clr'),
  groups = cget(groups,'clr'),
  steps  = cget(steps,'clr'),
  case   = cget(cases,'clr'))
lab = list(
  pop    = cget(pops,'lab'),
  part   = cget(parts,'lab'),
  part   = cget(parts,'lab'),
  groups = cget(groups,'lab'),
  steps  = cget(steps,'lab'),
  case   = cget(cases,'lab'))
csv.name = function(phase,key,case,b){
  return(root.path('data','mid',uid,sprintf('%d',N$cal),paste0(phase,'_',key,'_',case,'_',b,'.csv')))
}
load.csvs = function(phase,key,cases.alt,b='all'){
  if (missing(cases.alt)){ cases.alt = cases }
  X = do.call(rbind,lapply(names(cases.alt),function(case){
    X.i = do.call(rbind,lapply(b,function(bi){
      return(read.csv(csv.name(phase,key,case,b=bi)))
    }))
    X.i$case = case
    return(X.i)
  }))
  X$case.lab = factor(X$case,levels=names(cases.alt),labels=cget(cases.alt,'lab'))
  X$case.id  = factor(X$case,levels=names(cases.alt),labels=cget(cases.alt,'id'))
  return(X)
}
