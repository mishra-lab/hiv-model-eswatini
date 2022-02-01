suppressPackageStartupMessages({
  # library('rms')
  library('reshape2')
  library('ggplot2')
  library('viridis')
})
source('utils/ops.r')
N    = 1024
Nr   = 10
uid  = '2022-01-24'
load.data = function(){
  return(do.call(rbind,lapply(c('base',paste0('RL',Nr)),function(case){
    return(read.csv(root.path('data','mid',uid,paste0('sens_',case,'_N=0-',N-1,'.csv'))))
  })))
}
clean.data = function(X.raw,t.hor=2040){
  X.base = X.raw[X.raw$case=='base',]
  X      = X.raw[X.raw$case!='base',]
  X$r.id = 1:Nr
  X = X[order(X$r.id,X$seed),]
  o = function(Xi,output){
    out = Xi[[paste0(output,'_',t.hor)]]
  }
  X$inc.red = (o(X,'incid')  - o(X.base,'incid'))  / o(X,'incid')
  X$inc.ext = (o(X,'incid')  - o(X.base,'incid'))  / o(X.base,'incid')
  X$inf.red = (o(X,'cuminf') - o(X.base,'cuminf')) / o(X,'cuminf')
  X$inf.ext = (o(X,'cuminf') - o(X.base,'cuminf')) / o(X.base,'cuminf')
  X$d.vls.all = X.base$vls_u_ALL - X$vls_u_ALL
  X$d.vls.fsw = X.base$vls_u_FSW - X$vls_u_FSW
  for (output in c('incid','cuminf')){ X[grepl(output,names(X))] = NULL } # remove raw outputs
  return(X)
}
do.plot = function(X,y='inf.red',ylab='% Infections Averted',...){
  X.long = melt(X,m=c('d.vls.all','d.vls.fsw'))
  levels(X.long$variable) = c('NVS Overall','NVS FSW')
  g = ggplot(X.long,aes_string(x='100 * value',y=paste('100 *',y),...)) +
    # geom_line(aes(group=seed),alpha=.3) +
    geom_point() +
    facet_grid(cols=vars(variable)) +
    scale_colour_viridis(option='inferno',begin=.1,end=.9) +
    labs(x='% Not Virally Suppressed (NVS)',y=ylab) +
    theme_light()
  return(g)
}
X = clean.data(load.data(),2020)
g = do.plot(X,color='prev_all')
ggsave('Rplots.pdf',w=6,h=3)
