# config stuff for ad-hoc analyses of model outputs

source('utils/ops.r')
source('utils/plot.r')

N = list(batch=100,hsam=1000,isam=100,imis=100) # h1000i100b100
# N = list(batch=50,hsam=1000,isam=100,imis=25) # h1000i25b50
# N = list(batch=5,hsam=100,isam=15,imis=15) # h100i15b5
uid = '2024-03-12'
nid = sprintf('h%di%db%d',N$hsam,N$imis,N$batch)
strats = list(
  # populations
  'all'   = list(rgb(.40,.00,.40),'Overall'),
  'aq'    = list(rgb(.60,.20,.60),'Lower Risk'),
  'w'     = list(rgb(.80,.20,.20),'Women'),
  'wl'    = list(rgb(1.0,.60,.60),'Women Low'),
  'wm'    = list(rgb(.80,.40,.40),'Women Med'),
  'fsw.l' = list(rgb(.60,.20,.20),'LR FSW'),
  'fsw.h' = list(rgb(.40,.00,.00),'HR FSW'),
  'fsw'   = list(rgb(.60,.00,.00),'FSW'),
  'm'     = list(rgb(.20,.20,.80),'Men'),
  'ml'    = list(rgb(.60,.60,1.0),'Men Low'),
  'mm'    = list(rgb(.40,.40,.80),'Men Med'),  
  'cli.l' = list(rgb(.20,.20,.60),'LR Clients'),
  'cli.h' = list(rgb(.00,.00,.40),'HR Clients'),
  'cli'   = list(rgb(.00,.00,.60),'Clients'),
  # partnerships
  'msp'   = list(rgb(.26,.04,.41),'Main / Spousal'),
  'cas'   = list(rgb(.58,.15,.40),'Casual'),
  'swo'   = list(rgb(.87,.32,.23),'One-Off Sex Work'),
  'swr'   = list(rgb(.99,.65,.04),'Repeat Sex Work'),
  'swx'   = list(rgb(.93,.48,.13),'Sex Work Overall'),
  # scenarios
  'foi-rd'   = list('#00CC99','Rate-Duration','62'),
  'foi-ry'   = list('#0099CC','Rate-1-Year','22'),
  'foi-py'   = list('#CC00CC','Proportion-1-Year','4212'),
  'fsw-cli+' = list('#FF0033','FSW','62'),
  'fsw+cli-' = list('#0066CC','Clients','22'),
  'fsw-cli-' = list('#990099','Both','3212'),
  'fsw+cli+' = list('#FF9900','Neither','6212'),
  'sens'  = list('#999999','Sensitivity'),
  'base'  = list('#999999','Base','solid'))
strats$aqf  = strats$aqt  = strats$aq
strats$fswf = strats$fswt = strats$fsw
strats$clif = strats$clit = strats$cli
# combinations of the above
sets = list(
  base    = c('base'),
  sens    = c('sens'),
  foi     = c('foi-rd','foi-ry','foi-py','base'),
  art     = c('fsw-cli+','fsw+cli-','fsw-cli-','fsw+cli+','base'),
  pop.all = c('wl','wm','fsw.l','fsw.h','ml','mm','cli.l','cli.h'),
  pop.cal = c('all','w','m','fsw'),
  pop.art = c('all','aq','fsw','cli'),
  ptr     = c('msp','cas','swo','swr'))
# collect colours, labels, linetypes for plotting
f = function(x,i){ ifelse(len(x) < i,NA,x[[i]]) }
strat.cols = sapply(names(strats),function(id){ f(strats[[id]],1) })
strat.labs = sapply(names(strats),function(id){ f(strats[[id]],2) })
strat.lts  = sapply(names(strats),function(id){ f(strats[[id]],3) })
set.cols = sapply(sets,function(ids){ unname(sapply(ids,function(id){ strat.cols[[id]] })) })
set.labs = sapply(sets,function(ids){ unname(sapply(ids,function(id){ strat.labs[[id]] })) })
set.lts  = sapply(sets,function(ids){ unname(sapply(ids,function(id){ strat.lts[[id]] })) })
set.labs$foi[4] = 'Effective Partners Adjustment'; set.cols$foi[4] = '#FF9900'

csv.name = function(phase,key,case,b='all',ext='.csv',log=''){
  # generate standardized .csv or .rdata filename
  fname = root.path('data','csv',uid,nid,paste0( phase,'_',key,'_',case,'_',b,ext))
  if (nchar(log)){ pout(log,': ',fname) }
  return(fname)
}
read.csvs = function(phase,key,set,b='all',skip=NULL,rdata=''){
  # read & rbind csv files, possibly save / load from rdata (much faster)
  if (rdata=='load'){ load(file=csv.name(phase,key,set,ext='.rdata',log='load')); return(X) }
  X = rbind.lapply(sets[[set]],function(case){
    if (case %in% skip){ return(NULL) }
    X.i = rbind.lapply(b,function(bi){
      read.csv(csv.name(phase,key,case,b=bi,log='load'),as.is=TRUE)
    })
    X.i$case = case
    return(X.i)
  })
  X$case.lab = factor(X$case,levels=sets[[set]],labels=set.labs[[set]])
  if (rdata=='save'){ save(X,file=csv.name(phase,key,set,ext='.rdata',log='save')) }
  return(X)
}

grep.i.col = function(X.wide){
  # select main data columns from out.expo with mode='id'
  grep('^i\\d+\\.\\d+\\.\\d+$',colnames(X.wide))
}

melt.expo.i = function(X.wide,...){
  # melt main data columns from out.expo with model='id'
  X = melt(filter.cols(X.wide,...),measure=grep.i.col(X.wide),var='id')
  return(X[order(X$case),])
}

qs = c(0,.025,.05,.1,.25,.4,.45,.475,.5,.525,.55,.6,.75,.9,.95,.975,1)
q3 = c(.025,.5,.975)

expo.qs = function(X,q=qs,trans=identity){
  # compute quantiles from melted out.expo data
  vars = colnames(X)[!grepl('^id$|^value$|^ss$',colnames(X))]
  f = formula(paste('value ~',paste(vars,collapse='+')))
  fun = function(x){ trans(quantile(x,p=q)) }
  X.q = do.call(data.frame,aggregate(f,X,fun))
  colnames(X.q)[grep('^value\\.',colnames(X.q))] = paste0('q',q)
  return(X.q)
}
