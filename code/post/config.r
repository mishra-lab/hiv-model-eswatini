source('utils/ops.r')
source('utils/plot.r')

uid = '2022-04-20'
N = list(sam=100000,batchs=10,sens=10,topfit=.01)
spec = list(
pops = list(
  'wl'    = list(clr=rgb(1.,.6,.6),lab='Women Low'),
  'wm'    = list(clr=rgb(.8,.4,.4),lab='Women Med'),
  'fsw.l' = list(clr=rgb(.6,.2,.2),lab='LR FSW'),
  'fsw.h' = list(clr=rgb(.4,.0,.0),lab='HR FSW'),
  'ml'    = list(clr=rgb(.6,.6,1.),lab='Men Low'),
  'mm'    = list(clr=rgb(.4,.4,.8),lab='Men Med'),
  'cli.l' = list(clr=rgb(.2,.2,.6),lab='LR Clients'),
  'cli.h' = list(clr=rgb(.0,.0,.4),lab='HR Clients')),
groups = list(
  'all' = list(clr=rgb(.4,.4,.4),lab='Overall'),
  'aq'  = list(clr=rgb(.9,.7,.0),lab='Lower Risk'),
  'fsw' = list(clr=rgb(.6,.0,.0),lab='FSW'),
  'cli' = list(clr=rgb(.0,.0,.6),lab='Clients')),
parts = list(
  'msp' = list(clr=rgb(.26,.04,.41),lab='Main / Spousal'),
  'cas' = list(clr=rgb(.58,.15,.40),lab='Casual'),
  'swo' = list(clr=rgb(.87,.32,.23),lab='Occas. Sex Work'),
  'swr' = list(clr=rgb(.99,.65,.04),lab='Regular Sex Work')),
t.hors = list(
  '2020' = list(clr=rgb(.26,.04,.41),lab='2020'),
  '2030' = list(clr=rgb(.69,.20,.35),lab='2030'),
  '2040' = list(clr=rgb(.98,.55,.04),lab='2040')),
cascade = list(
  'diagnosed' = list(clr=rgb(.0,.0,.0),lab='Diagnosed among HIV+'),
  'treated.c' = list(clr=rgb(.0,.0,.0),lab='Treated among diagnosed'),
  'vls.c'     = list(clr=rgb(.0,.0,.0),lab='VLS among treated'),
  'treated.u' = list(clr=rgb(.0,.0,.0),lab='Treated among HIV+'),
  'vls.u'     = list(clr=rgb(.0,.0,.0),lab='VLS among HIV+')),
cases.art = list(
  'fsw-cli+' = list(clr=rgb(.9,.3,.3),id='-+',lab='Left Behind: FSW'),
  'fsw+cli-' = list(clr=rgb(.3,.3,.9),id='+-',lab='Left Behind: Clients'),
  'fsw-cli-' = list(clr=rgb(.8,.3,.8),id='--',lab='Left Behind: FSW & Clients'),
  'fsw+cli+' = list(clr=rgb(.9,.7,.0),id='++',lab='Left Behind: Neither'),
  'base'     = list(clr=rgb(.4,.4,.4),id='bc',lab='Base Case')),
cases.foi = list(
  # 'lin'  = list(clr=rgb(.500,.500,.500),id='lin',lab='<3>'),
  'bpd'  = list(clr=rgb(.267,.005,.329),id='bpd',lab='<1a>'),
  'bpy'  = list(clr=rgb(.191,.407,.556),id='bpy',lab='<1b>'),
  'bmy'  = list(clr=rgb(.208,.719,.473),id='bmy',lab='<2b>'),
  'base' = list(clr=rgb(.993,.906,.144),id='base',lab='<4*>')))
sget = function(name,subname){
  return(unname(sapply(spec[[name]],function(def){def[[subname]]})))
}
gen.name = function(phase,key,case,b='all',ext='.csv'){
  return(root.path('data','csv',uid,sprintf('%d',N$sam),
    paste0( phase,'_',key,'_',case,'_',b,ext)))
}
read.csvs = function(phase,key,name,b='all',ext='.csv',skip=NULL){
  cases = spec[[name]]
  X = do.call(rbind,lapply(names(cases),function(case){
    if (case %in% skip){ return(NULL) }
    X.i = do.call(rbind,lapply(b,function(bi){
      return(read.csv(gen.name(phase,key,case,b=bi,ext=ext)))
    }))
    X.i$case = case
    return(X.i)
  }))
  X$case.lab = factor(X$case,levels=names(cases),labels=sget(name,'lab'))
  X$case.id  = factor(X$case,levels=names(cases),labels=sget(name,'id'))
  return(X)
}