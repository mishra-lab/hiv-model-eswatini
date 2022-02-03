suppressPackageStartupMessages({
  library('geepack')
  library('reshape2')
  library('ggplot2')
  library('viridis')
})
source('utils/ops.r')
# config: model fits
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
  X$ipr_all   = X$inc_all / X$prev_all
  for (output in c('incid','cuminf')){ X[grepl(output,names(X))] = NULL } # remove raw outputs
  return(X)
}
do.glm = function(X,y='inf.red',pred.vars,mod.vars,scale=FALSE){
  f = paste(y,' ~ ',paste(sep=' + ',0,
    paste(pred.vars,collapse=' * '),
    paste('(',paste(pred.vars,collapse=' + '),') : (',paste(mod.vars,collapse=' + '),')')))
  if (scale){
    vars = c(pred.vars,mod.vars)
    X[,vars] = apply(X[,vars],2,function(x){ (x-mean(x))/sd(x) })
  }
  print(f)
  m = geeglm(formula(f),'gaussian',X,id=factor(X$seed),corstr='i')
  return(m)
}
plot.effects = function(models,c.lab='Year'){
  w = qnorm(.975)
  E = do.call(rbind,lapply(names(models),function(name){
    e = summary(models[[name]])$coefficients
    e$vars = rownames(e)
    e$est.low  = e$Estimate - w*e$Std.err
    e$est.high = e$Estimate + w*e$Std.err
    e$name = name
    return(e)
  }))
  g = ggplot(E,aes(y=vars,x=Estimate,xmin=est.low,xmax=est.high,color=name)) +
    geom_point(position=position_dodge(width=.8)) +
    geom_linerange(,position=position_dodge(width=.8)) +
    scale_colour_viridis_d(option='inferno',begin=.2,end=.8) +
    scale_x_continuous(trans='pseudo_log') +
    labs(x='Effect',y='Condition',color=c.lab) +
    theme_light()
  return(g)
}
plot.points = function(X,y,ylab,...){
  # TODO: facet by effects?
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
fig.save = function(...,w=7,h=7){
  fig.name = root.path('out','fig',uid,paste0(paste(...,sep='_'),'.pdf'))
  print(paste('saving:',fig.name))
  ggsave(fig.name,w=w,h=h)
}
# config: this analysis
y.def = 'inf.red'
y.lab = '% Infections Averted'
t.def = 2040
t.vec = c(2010,2020,2030,2040,2050)
pred.vars = c('d.vls.all','d.vls.fsw')
mod.vars = list(
  mech = c('EHY_acute','RPA_condom_a.v','PX_fsw','dur_fsw_l','A_swq_cli','dur_cli'),
  cond = c('prev_ratio_fsw.wq','prev_ratio_cli.mq','inc_all','prev_all')
)
X.raw = load.data()
X = clean.data(X.raw,t.def)
for (name in names(mod.vars)){
  # standardized effects
  M = lapply(setNames(t.vec,t.vec),function(t){
    do.glm(clean.data(X.raw,t),y.def,pred.vars,mod.vars[[name]],scale=TRUE)
  })
  g = plot.effects(M)
  fig.save('sens','std-eff',name,w=6,h=8)
  # raw effects
  m = do.glm(X,y.def,pred.vars,mod.vars[[name]],scale=FALSE)
  print(summary(m))
}
# exploratory
# g = plot.points(X,y.def,y.lab,color='EHY_acute')
# ggsave('Rplots.pdf',w=8,h=4)
# TODO: clean up labels...