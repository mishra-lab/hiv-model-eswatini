library('ggplot2')
library('reshape2')
source('utils/ops.r')
source('utils/plot.r')
# --------------------------------------------------------------------------------------------------
# config: model fits
N    = 1024
uid  = '2022-02-10'
# config: cases
cases = list(
  'FSW-Cli+AQ+'=list(clr=rgb(.9,.3,.3),cid='-++',title='Leave Behind: FSW'),
  'FSW+Cli-AQ+'=list(clr=rgb(.3,.3,.9),cid='+-+',title='Leave Behind: Clients'),
  'FSW+Cli+AQ-'=list(clr=rgb(.9,.9,.3),cid='++-',title='Leave Behind: Lower Risk'),
  'FSW-Cli-AQ+'=list(clr=rgb(.8,.3,.8),cid='--+',title='Leave Behind: FSW & Clients'),
  'FSW-Cli+AQ-'=list(clr=rgb(.9,.6,.3),cid='-+-',title='Leave Behind: FSW & Lower Risk'),
  'FSW+Cli-AQ-'=list(clr=rgb(.3,.9,.6),cid='+--',title='Leave Behind: Clients & Lower Risk'),
  'FSW-Cli-AQ-'=list(clr=rgb(.5,.5,.5),cid='---',title='Worst Case'),
  'FSW+Cli+AQ+'=list(clr=rgb(.8,.8,.8),cid='+++',title='Best Case'))
colors = unname(sapply(cases,function(case){case$clr}))
titles = unname(sapply(cases,function(case){case$title}))
cids   = unname(sapply(cases,function(case){case$cid}))
t.hors = seq(2000,2050,5)
t.labs = t.hors; t.labs[seq(2,length(t.labs),2)] = '' # HACK
# --------------------------------------------------------------------------------------------------
# load & process data
load.outs.data = function(){
  X = do.call(rbind,lapply(names(cases),function(case){
    X.i = read.csv(root.path('data','mid',uid,paste0('outs_',case,'_N=0-',N-1,'.csv')))
    return(X.i)
  }))
  X$title = X$case; levels(X$title) = titles
  X$cid   = X$case; levels(X$cid)   = cids
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
  X.lo = X[X$cid=='---',]
  X.hi = X[X$cid=='+++',]
  # X$PX_cli  = X$PX_cli / .48 # HACK TODO
  X$inf.red = (X.lo$cuminf - X$cuminf) / X.lo$cuminf
  X$inf.ext = (X$cuminf - X.hi$cuminf) / X.hi$cuminf
  X$inc.red = (X.lo$incid  - X$incid)  / X.lo$incid
  X$inc.ext = (X$incid  - X.hi$incid)  / X.hi$incid
  return(X)
}
plot.out = function(X,ylist,ylab,scale=1,tlims=NULL,ncol=NULL,...){
  args = list(...)
  X.long = melt(X,m=unlist(ylist))
  levels(X.long$variable) = names(ylist)
  if (!is.null(tlims)){ t = X.long$t; X.long = X.long[t>=tlims[1] & t<=tlims[2],] }
  g = ggplot(X.long,aes_string(x='t',y=paste(scale,'*value'),...)) +
    geom_line() +
    labs(x='Year',y=ylab) +
    theme_light()
  if (length(ylist) > 1){ g = g + facet_wrap(vars(variable),ncol=ncol,scales='free_y') }
  if ('color' %in% names(args)){ g = g + scale_color_manual(values=colors) }
  if ('fill' %in% names(args)){ g = g + scale_fill_manual(values=colors) }
  return(g)
}
plot.obj.1 = function(X,y='100*inf.red',ylab='Infections averted (%)',...){
  if (grepl('\\.red',y)){ skip = '---'; ord = match(c('--+','-++','+-+','+++','-+-','+--','++-','---'),cids) }
  if (grepl('\\.ext',y)){ skip = '+++'; ord = match(c('--+','-++','+-+','---','-+-','+--','++-','+++'),cids) }
  X$retitle = factor(X$title,levels = titles[ord])
  g = ggplot(X[X$cid!=skip,],aes_string(x='factor(year)',y=y,fill='retitle',...)) +
    geom_boxplot(show.legend=FALSE,outlier.size=.5,lwd=.3,width=.6) +
    facet_wrap(vars(retitle),ncol=3,dir='v') +
    labs(x='Year',y=ylab) +
    scale_x_discrete(labels=t.labs) +
    scale_fill_manual(values=colors[ord]) +
    scale_color_manual(values=colors[ord]) +
    theme_light()
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
  Nf = length(unique(X.fun('+++','seed')))
  save.tex(N,'n.sample',dec=0)
  save.tex(10,'top.pct.fit',dec=0)
  save.tex(Nf,'n.fit',dec=0)
  save.tex(10,'n.rand',dec=0) # obj.2
  save.med.ci(X.fun('+++','PX_fsw'),'px.fsw',per=100)
  save.med.ci(X.fun('+++','PX_cli'),'px.cli',per=100)
  for (cid in cids){
    save.med.ci(X.fun(cid,'prev_all'),         paste0(cid,'/prev.all.2020'),per=100)
    save.med.ci(X.fun(cid,'inc_all'),          paste0(cid,'/inc.all.2020'),per=1000)
    save.med.ci(X.fun(cid,'prev_ratio_fsw.wq'),paste0(cid,'/pr.fsw.2020'),dec=2)
    save.med.ci(X.fun(cid,'prev_ratio_cli.mq'),paste0(cid,'/pr.cli.2020'),dec=2)
    for (year in seq(2000,2050,10)){
      save.med.ci(X.fun(cid,'inf.red',year=year),paste0(cid,'/inf.red.',year),per=100)
      save.med.ci(1 - (X.fun(cid,'inf.red',year=year)/X.fun('+++','inf.red',year=year)),
        paste0(cid,'/inf.red.vs.best.',year),per=100)
    }
  }
}
# --------------------------------------------------------------------------------------------------
Y = load.outs.data()
ylist = list(
  'Overall'    = 'inc_all.mu',
  'Lower Risk' = 'inc_aq.mu',
  'FSW'        = 'inc_fsw.mu',
  'Clients'    = 'inc_cli.mu')
g = plot.out(Y,ylist,ylab='HIV Incidence (per 1000 PY)',scale=1000,color='title',tlims=c(2000,2050)) + labs(color='')
ggsave('Rplots.pdf',w=8,h=4)
X = load.sens.data()
numeric.obj.1(X);
plot.obj.1(X,y='100*inf.red',ylab='Infections averted (%) vs worst case');      fig.save(uid,'obj_1_inf_red')
plot.obj.1(X,y='100*inf.ext',ylab='Additional infections (%) vs best case');    fig.save(uid,'obj_1_inf_ext')
plot.obj.1(X,y='100*inc.red',ylab='Incidence reduction (%) vs worst case');     fig.save(uid,'obj_1_inc_red')
plot.obj.1(X,y='  1*inc.ext',ylab='Additional incidence (times) vs best case'); fig.save(uid,'obj_1_inc_ext')

