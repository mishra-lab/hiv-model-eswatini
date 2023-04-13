source('utils/ops.r')
source('utils/plot.r')

uid = '2023-03-03'
N = list(sam=100000,batchs=10,sens=10,topfit=.01)
slicers = list(
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
  'msp'   = list(rgb(.26,.04,.41),'Main / Spousal'),
  'cas'   = list(rgb(.58,.15,.40),'Casual'),
  'swo'   = list(rgb(.87,.32,.23),'Sex Work Occas.'),
  'swr'   = list(rgb(.99,.65,.04),'Sex Work Reg.'),
  'swx'   = list(rgb(.93,.48,.13),'Sex Work Overall'),
  'rd'    = list('#FF0033','Rate-Duration','62'),
  'ry'    = list('#990099','Rate-1-Year','22'),
  'pd'    = list('#00CC66','Proportion-Duration','8212'),
  'py'    = list('#0066CC','Proportion-1-Year','4212'),
  'fsw-cli+' = list('#FF0033','FSW'),
  'fsw+cli-' = list('#0066CC','Clients'),
  'fsw-cli-' = list('#990099','Both'),
  'fsw+cli+' = list('#FF9900','Neither'),
  'base'  = list('#999999','Base','solid'))
slicers$aqf  = slicers$aqt  = slicers$aq
slicers$fswf = slicers$fswt = slicers$fsw
slicers$clif = slicers$clit = slicers$cli
sets = list(
  base    = c('base'),
  foi     = c('rd','ry','pd','py','base'),
  art     = c('fsw-cli+','fsw+cli-','fsw-cli-','fsw+cli+','base'),
  pop.all = c('wl','wm','fsw.l','fsw.h','ml','mm','cli.l','cli.h'),
  pop.cal = c('all','w','m','fsw'),
  pop.art = c('all','aq','fsw','cli'),
  part    = c('msp','cas','swo','swr'))
slice.cols = sapply(names(slicers),function(id){ slicers[[id]][[1]] })
slice.labs = sapply(names(slicers),function(id){ slicers[[id]][[2]] })
set.cols = sapply(sets,function(ids){ unname(sapply(ids,function(id){ slice.cols[[id]] })) })
set.labs = sapply(sets,function(ids){ unname(sapply(ids,function(id){ slice.labs[[id]] })) })
set.lts = list(foi=unname(sapply(sets$foi,function(id){ slicers[[id]][[3]] }))) # TODO: better way?
set.labs$foi[5] = 'New Proposed'; set.cols$foi[5] = '#FF9900'

gen.name = function(phase,key,case,b='all',ext='.csv'){
  return(root.path('data','csv',uid,sprintf('%d',N$sam),
    paste0( phase,'_',key,'_',case,'_',b,ext)))
}
read.csvs = function(phase,key,set,b='all',skip=NULL,rdata=''){
  if (rdata=='load'){ load(file=gen.name(phase,key,set,ext='.rdata')); return(X) }
  X = do.call(rbind,lapply(sets[[set]],function(case){
    if (case %in% skip){ return(NULL) }
    X.i = do.call(rbind,lapply(b,function(bi){
      read.csv(gen.name(phase,key,case,b=bi),as.is=TRUE)
    }))
    X.i$case = case
    return(X.i)
  }))
  X$case.lab = factor(X$case,levels=sets[[set]],labels=set.labs[[set]])
  if (rdata=='save'){ save(X,file=gen.name(phase,key,set,ext='.rdata')) }
  return(X)
}

melt.expo.s = function(X,...){
  X = melt(filter.cols(X,...),measure=grep('^s\\d*',colnames(X)),var='seed')
  X$seed = as.integer(gsub('s','',X$seed))
  return(X[order(X$case),])
}

qs = c(0,.025,.05,.1,.25,.4,.45,.475,.5,.525,.55,.6,.75,.9,.95,.975,1)

expo.qs = function(X,q=qs){
  vars = colnames(X)[!grepl('^seed$|^value$',colnames(X))]
  f = formula(paste('value ~',paste(vars,collapse='+')))
  X.q = do.call(data.frame,aggregate(f,X,quantile,p=q))
  colnames(X.q)[grep('^value\\.',colnames(X.q))] = paste0('q',q)
  return(X.q)
}
