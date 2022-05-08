source('post/config.r')

# --------------------------------------------------------------------------------------------------
# load & process data
load.keyout.data = function(t.hors){
  if (missing(t.hors)){ t.hors = seq(2000,2050,5) }
  X.raw = read.csvs('art','keyout','cases.art')
  vars = c('cuminf_all','inc_all','prev_all')
  X = reshape(X.raw,idvar=c('seed','case'),direction='long',times=t.hors,timevar='year',
    v.names=vars,varying=lapply(vars,function(v){ grepl(v,names(X.raw)) }))
  X = X[order(X$case,X$seed),]
  X.base = X[X$case.id=='bc',]
  X$inf.red  = (X$cuminf_all - X.base$cuminf_all) / X$cuminf_all
  X$inf.add  = (X$cuminf_all - X.base$cuminf_all) / X.base$cuminf_all
  X$inc.red  = (X$inc_all    - X.base$inc_all)    / X$inc_all
  X$inc.add  = (X$inc_all    - X.base$inc_all)    / X.base$inc_all
  X$prev.red = (X$prev_all   - X.base$prev_all)   / X$prev_all
  X$prev.add = (X$prev_all   - X.base$prev_all)   / X.base$prev_all
  return(X)
}
plot.out = function(X,ylist,ylab,scale=1,tlims=NULL,...){
  args = list(...)
  X = aggregate(.~case+t+case.lab+case.id,X,median) # if >1 batches
  X.long = melt(X,m=unlist(ylist))
  levels(X.long$variable) = names(ylist)
  if (!is.null(tlims)){ t = X.long$t; X.long = X.long[t>=tlims[1] & t<=tlims[2],] }
  g = ggplot(X.long,aes_string(x='t',y=paste(scale,'*value'),...)) +
    geom_line() +
    labs(x='Year',y=ylab) +
    theme_light()
  if (length(ylist) > 1){ g = g + facet_wrap(vars(variable),nrow=1,scales='free_y') }
  if ('color' %in% names(args)){ g = g + scale_color_manual(values=sget('cases.art','clr')) }
  if ('fill' %in% names(args)){ g = g + scale_fill_manual(values=sget('cases.art','clr')) }
  return(g)
}
plot.obj.1 = function(X,y='100*inf.red',ylab='Infections averted (%)',yl,tlims=c(2005,2040),...){
  X$case.lab = factor(X$case.lab,labels=gsub('Left Behind: ','',levels(X$case.lab)))
  if (!is.null(tlims)){ t = as.numeric(X$year); X = X[t>=tlims[1] & t<=tlims[2],] }
  g = ggplot(X[X$case.id!='bc',],aes_string(x='factor(year)',y=y,...)) +
    geom_boxplot(aes(color=case.lab),outlier.size=.5,lwd=.0,width=.6,position=position_dodge(.8)) +
    geom_boxplot(aes(fill=case.lab),outlier.color=NA,lwd=.3,width=.6,position=position_dodge(.8)) +
    labs(x='Year',y=ylab,color='Left Behind:',fill='Left Behind:') + lims(y=yl) +
    scale_fill_manual(values=sget('cases.art','clr')) +
    scale_color_manual(values=sget('cases.art','clr')) +
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
  v = quantile(values,c(.025,.5,.975),na.rm=na.rm)
  save.tex(v[2],paste0(fname,'.md'),dec=dec,per=per)
  save.tex(v[1],paste0(fname,'.lo'),dec=dec,per=per)
  save.tex(v[3],paste0(fname,'.hi'),dec=dec,per=per)
}
numeric.obj.1 = function(X){
  X.fun = function(case.id,col,year=NULL){
    if (is.null(year)){ Xi = X } else { Xi = X[X$year==year,] }
    return(Xi[Xi$case.id==case.id,][[col]])
  }
  X$prev_ratio_fsw.wq_2020 = X$prev_fsw_2020 / X$prev_wq_2020
  X$prev_ratio_cli.mq_2020 = X$prev_cli_2020 / X$prev_mq_2020
  save.tex(N$sam,'n.sam',dec=0)
  save.tex(100*N$topfit,'top.pct.fit',dec=0)
  save.tex(N$sam*N$topfit,'n.fit',dec=0)
  save.tex(N$sens,'n.sens.per',dec=0) # obj.2
  save.tex(N$sens*N$topfit,'n.sens',dec=0) # obj.2
  save.med.ci(X.fun('bc','PX_fsw'),'px.fsw',per=100)
  save.med.ci(X.fun('bc','PX_cli'),'px.cli',per=100)
  save.med.ci(X.fun('bc','Rdx_fsw'),'Rdx.fsw',dec=2)
  save.med.ci(X.fun('bc','Rdx_cli'),'Rdx.cli',dec=2)
  save.med.ci(X.fun('bc','Rtx_fsw'),'Rtx.fsw',dec=2)
  save.med.ci(X.fun('bc','Rtx_cli'),'Rtx.cli',dec=2)
  for (cid in sget('cases.art','id')){
    save.med.ci(X.fun(cid,'prev_all',year=2020),    paste0(cid,'/prev.all.2020'),per=100)
    save.med.ci(X.fun(cid,'inc_all',year=2020),     paste0(cid,'/inc.all.2020'),per=1000)
    save.med.ci(X.fun(cid,'prev_ratio_fsw.wq_2020'),paste0(cid,'/pr.fsw.2020'),dec=2)
    save.med.ci(X.fun(cid,'prev_ratio_cli.mq_2020'),paste0(cid,'/pr.cli.2020'),dec=2)
    for (year in seq(2000,2050,10)){
      save.med.ci(X.fun(cid,'inf.add',year=year),paste0(cid,'/inf.add.',year),per=100)
      save.med.ci(1 - (X.fun(cid,'inf.add',year=year)/X.fun('--','inf.add',year=year)),
        paste0(cid,'/inf.add.vs.--.',year),per=100)
    }
  }
}
# --------------------------------------------------------------------------------------------------
# TODO: do we really need this ... ? maybe just plot in python
# X.out = load.csvs('outs')
# ylist = list('Overall'='inc_all.mu','Lower Risk'='inc_aq.mu','FSW'='inc_fsw.mu','Clients'='inc_cli.mu')
# g = plot.out(X.out,ylist,ylab='HIV Incidence (per 1000 PY)',scale=1000,color='case.lab',tlims=c(2000,2050)) +
#   labs(color='Left Behind:') + theme(legend.position='top'); fig.save(uid,'obj_1_inc',w=10,h=3)
X = load.keyout.data()
numeric.obj.1(X)
plot.obj.1(X,y='100*inf.red',ylab='Infections averted (%)',      yl=c(0, 60)); fig.save(uid,'art_1_inf_red',w=8,h=4)
plot.obj.1(X,y='100*inf.add',ylab='Additional infections (%)',   yl=c(0,130)); fig.save(uid,'art_1_inf_add',w=8,h=4)
plot.obj.1(X,y='100*inc.red',ylab='Incidence reduction (%)',     yl=c(0,100)); fig.save(uid,'art_1_inc_red',w=8,h=4)
plot.obj.1(X,y='  1*inc.add',ylab='Additional incidence (times)',yl=c(0, 30)); fig.save(uid,'art_1_inc_add',w=8,h=4)

