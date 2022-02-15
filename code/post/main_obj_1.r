library('ggplot2')
library('reshape2')
source('utils/ops.r')
source('utils/plot.r')
# --------------------------------------------------------------------------------------------------
# config: model fits
N    = 1024
uid  = '2022-02-13'
# config: cases
cases = list(
  'FSW-Cli+'=list(clr=rgb(.9,.3,.3),cid='-+',title='FSW'),
  'FSW+Cli-'=list(clr=rgb(.3,.3,.9),cid='+-',title='Clients'),
  'FSW-Cli-'=list(clr=rgb(.8,.3,.8),cid='--',title='FSW & Clients'),
  'FSW+Cli+'=list(clr=rgb(.9,.7,.0),cid='++',title='Neither'),
  'base'    =list(clr=rgb(.4,.4,.4),cid='bc',title='Base Case'))
colors = unname(sapply(cases,function(case){case$clr}))
titles = unname(sapply(cases,function(case){case$title}))
cids   = unname(sapply(cases,function(case){case$cid}))
t.hors = seq(2000,2050,5)
# --------------------------------------------------------------------------------------------------
# load & process data
load.outs.data = function(){
  X = do.call(rbind,lapply(names(cases),function(case){
    X.i = read.csv(root.path('data','mid',uid,paste0('outs_',case,'_N=0-',N-1,'.csv')))
    return(X.i)
  }))
  X$title = factor(X$case,levels=names(cases),labels=titles)
  X$cid   = factor(X$case,levels=names(cases),labels=cids)
  return(X)
}
load.sens.data = function(){
  X.raw = do.call(rbind,lapply(names(cases),function(case){
    X.raw.i = read.csv(root.path('data','mid',uid,paste0('sens_',case,'_N=0-',N-1,'.csv')))
    X.raw.i$label = cases[[case]]$label
    return(X.raw.i)
  }))
  vars = c('cuminf','incid')
  X = reshape(X.raw,idvar=c('seed','case'),direction='long',times=t.hors,timevar='year',
    v.names=vars,varying=lapply(vars,function(i){ grepl(i,names(X.raw)) }))
  X = X[order(X$case,X$seed),]
  X$title = X$case; levels(X$title) = titles
  X$cids  = X$case; levels(X$cid)   = cids
  X.base = X[X$cid=='bc',]
  X$inf.red = (X$cuminf - X.base$cuminf) / X$cuminf
  X$inf.add = (X$cuminf - X.base$cuminf) / X.base$cuminf
  X$inc.red = (X$incid  - X.base$incid)  / X$incid
  X$inc.add = (X$incid  - X.base$incid)  / X.base$incid
  return(X)
}
plot.out = function(X,ylist,ylab,scale=1,tlims=NULL,...){
  args = list(...)
  X.long = melt(X,m=unlist(ylist))
  levels(X.long$variable) = names(ylist)
  if (!is.null(tlims)){ t = X.long$t; X.long = X.long[t>=tlims[1] & t<=tlims[2],] }
  g = ggplot(X.long,aes_string(x='t',y=paste(scale,'*value'),...)) +
    geom_line() +
    labs(x='Year',y=ylab) +
    theme_light()
  if (length(ylist) > 1){ g = g + facet_wrap(vars(variable),nrow=1,scales='free_y') }
  if ('color' %in% names(args)){ g = g + scale_color_manual(values=colors) }
  if ('fill' %in% names(args)){ g = g + scale_fill_manual(values=colors) }
  return(g)
}
plot.obj.1 = function(X,y='100*inf.red',ylab='Infections averted (%)',yl,...){
  g = ggplot(X[X$cid!='bc',],aes_string(x='factor(year)',y=y,...)) +
    geom_boxplot(aes(color=title),outlier.size=.5,lwd=.0,width=.6,position=position_dodge(.8)) +
    geom_boxplot(aes(fill=title),outlier.color=NA,lwd=.3,width=.6,position=position_dodge(.8)) +
    labs(x='Year',y=ylab,color='Left Behind:',fill='Left Behind:') + lims(y=yl) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    theme_light() +
    theme(legend.position=c(.01,.99),legend.justification=c(0,1))
  return(g)
}
save.tex = function(value,...,dec=3,per=1){
  value = value * per
  dec = dec - log10(per)
  str = sprintf(paste0('%.',dec,'f'),value)
  writeLines(str,root.path('out','tex',uid,...,create=TRUE))
}
save.med.ci = function(values,fname,dec=3,per=1,na.rm=TRUE){
  v = quantile(values,c(.025, .5, .975),na.rm=na.rm)
  save.tex(v[2],paste0(fname,'.md'),dec=dec,per=per)
  save.tex(v[1],paste0(fname,'.lo'),dec=dec,per=per)
  save.tex(v[3],paste0(fname,'.hi'),dec=dec,per=per)
}
numeric.obj.1 = function(X){
  X.fun = function(cid,col,year=NULL){
    if (is.null(year)){ Xi = X } else { Xi = X[X$year==year,] }
    return(Xi[Xi$cid==cid,][[col]])
  }
  Nf = length(unique(X.fun('bc','seed')))
  save.tex(N,'n.sample',dec=0)
  save.tex(10,'top.pct.fit',dec=0)
  save.tex(Nf,'n.fit',dec=0)
  save.tex(10,'n.rand',dec=0) # obj.2
  save.med.ci(X.fun('bc','PX_fsw'),'px.fsw',per=100)
  save.med.ci(X.fun('bc','PX_cli'),'px.cli',per=100)
  save.med.ci(X.fun('bc','Rdx_fsw'),'Rdx.fsw',dec=2)
  save.med.ci(X.fun('bc','Rdx_cli'),'Rdx.cli',dec=2)
  save.med.ci(X.fun('bc','Rtx_fsw'),'Rtx.fsw',dec=2)
  save.med.ci(X.fun('bc','Rtx_cli'),'Rtx.cli',dec=2)
  for (cid in cids){
    save.med.ci(X.fun(cid,'prev_all'),         paste0(cid,'/prev.all.2020'),per=100)
    save.med.ci(X.fun(cid,'inc_all'),          paste0(cid,'/inc.all.2020'),per=1000)
    save.med.ci(X.fun(cid,'prev_ratio_fsw.wq'),paste0(cid,'/pr.fsw.2020'),dec=2)
    save.med.ci(X.fun(cid,'prev_ratio_cli.mq'),paste0(cid,'/pr.cli.2020'),dec=2)
    for (year in seq(2000,2050,10)){
      save.med.ci(X.fun(cid,'inf.add',year=year),paste0(cid,'/inf.add.',year),per=100)
      save.med.ci(1 - (X.fun(cid,'inf.add',year=year)/X.fun('--','inf.add',year=year)),
        paste0(cid,'/inf.add.vs.--.',year),per=100)
    }
  }
}
# --------------------------------------------------------------------------------------------------
Y = load.outs.data()
ylist = list('Overall'='inc_all.mu','Lower Risk'='inc_aq.mu','FSW'='inc_fsw.mu','Clients'='inc_cli.mu')
g = plot.out(Y,ylist,ylab='HIV Incidence (per 1000 PY)',scale=1000,color='title',tlims=c(2000,2050)) +
  labs(color='Left Behind:') + theme(legend.position='top'); fig.save(uid,'obj_1_inc',w=10,h=3)
X = load.sens.data()
numeric.obj.1(X)
plot.obj.1(X,y='100*inf.red',ylab='Infections averted (%)',      yl=c(0, 70)); fig.save(uid,'obj_1_inf_red',w=8,h=4)
plot.obj.1(X,y='100*inf.add',ylab='Additional infections (%)',   yl=c(0,150)); fig.save(uid,'obj_1_inf_add',w=8,h=4)
plot.obj.1(X,y='100*inc.red',ylab='Incidence reduction (%)',     yl=c(0,100)); fig.save(uid,'obj_1_inc_red',w=8,h=4)
plot.obj.1(X,y='  1*inc.add',ylab='Additional incidence (times)',yl=c(0, 30)); fig.save(uid,'obj_1_inc_add',w=8,h=4)

