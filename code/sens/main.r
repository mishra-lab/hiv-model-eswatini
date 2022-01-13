suppressPackageStartupMessages({
  # library('rms')
  library('epiR')
  library('ggplot2')
})
source('utils/ops.r')
uid  = '2022-01-10'
cases = c(
  'Low Overall\nLow FSW' ='LoLo',
  'Low Overall\nHigh FSW'='LoHi',
  'High Overall\nLow FSW'='HiLo')

load.data = function(){
  return(do.call(rbind,lapply(c('base',cases),function(case){
    return(read.csv(root.path('data','mid',uid,paste0('sens_',case,'_N=0-699.csv'))))
  })))
}
clean.data = function(X.raw,t.hor){
  X.base = X.raw[X.raw$case=='base',]
  X = X.raw[X.raw$case!='base',]
  o = function(output){ paste0(output,'_',t.hor) }
  X$inc.red = (X[[o('incid')]]  - X.base[[o('incid')]])  / X[[o('incid')]]
  X$inc.ext = (X[[o('incid')]]  - X.base[[o('incid')]])  / X.base[[o('incid')]]
  X$inf.red = (X[[o('cuminf')]] - X.base[[o('cuminf')]]) / X[[o('cuminf')]]
  X$inf.ext = (X[[o('cuminf')]] - X.base[[o('cuminf')]]) / X.base[[o('cuminf')]]
  for (output in c('incid','cuminf')){ X[grepl(output,names(X))] = NULL } # remove raw outputs
  return(X)
}
do.glm = function(X,y,vars){
  f = formula(paste(y,'~',paste(vars,collapse=' + ')))
  for (case in cases){
    M = Glm(f,X[X$case==case,],family=quasibinomial())
    print(M)
  }
}
do.prcc = function(X,y,vars){
  Q = do.call(rbind,lapply(cases,function(case){
    Q = epi.prcc(X[X$case==case,c(vars,y)]) # NOTE: order robust?
    Q$var  = vars
    Q$case = case
    Q$case.name = names(cases)[grepl(case,cases)]
    return(Q)
  }))
}
plot.prcc = function(Q,y){
  w = .7
  g.lab = 'Cascade Scenario\nvs Base Case:\nHigh Overall\nHigh FSW'
  g = ggplot(Q,aes(y=est,x=var,fill=case.name,color=case.name)) +
    geom_bar(stat='identity',width=w,position=position_dodge(w),alpha=.5) +
    geom_errorbar(aes(ymin=lower,ymax=upper),width=.2,position=position_dodge(w)) +
    labs(x='Parameter',y=paste0('PRCC (',y,')'),color=g.lab,fill=g.lab) +
    theme_light() +
    scale_fill_viridis_d(option='inferno',begin=.2,end=.8) +
    scale_colour_viridis_d(option='inferno',begin=.2,end=.8) +
    theme(legend.position='top',legend.justification='right') +
    coord_flip()
  return(g)
}

X.raw = load.data()
X = clean.data(X.raw,2050)
vars = setdiff(colnames(X),c('case','seed','inc.red','inf.red','inc.ext','inf.ext'))
# do.glm(X,'inf.ext',vars)
Q = do.prcc(X,'inf.ext',vars)
g = plot.prcc(Q,y='Relative Additional Infections')
ggsave('Rplots.pdf',h=2+.15*length(cases)*length(vars),w=5)