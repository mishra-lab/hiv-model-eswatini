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
load.expo.data = function(){
  E = read.csvs('art','expo','cases.art')
  E = E[E$pop %in% names(spec$groups),]
  E$out = gsub('_c','.c',gsub('_u','.u',E$out))
  E$step = factor(E$out,levels=names(spec$cascade),labels=sget('cascade','lab'))
  E$groups = factor(E$pop,levels=names(spec$groups),labels=sget('groups','lab'))
  return(E)
}
# --------------------------------------------------------------------------------------------------
# plot stuff
plot.out = function(E,outs,color=NULL,scale=1){
  g = ggplot(E[E$out %in% outs,],aes_string(x='t',color=color,fill=color)) +
    geom_ribbon(aes(ymin=scale*q0.05,ymax=scale*q0.95),alpha=.2,color=NA) +
    geom_line(aes(y=scale*q0.5)) +
    scale_color_manual(values=sget('cases.art','clr')) +
    scale_fill_manual(values=sget('cases.art','clr')) +
    theme_light()
}
plot.cascade = function(E,steps,y.lab='Cascade step (%)'){
  g = plot.out(E,steps,color='case.lab',scale=100) +
    facet_grid('step ~ groups') +
    labs(x='Year',y=y.lab,color='',fill='') +
    scale_x_continuous(limits=c(2000,2040),breaks=c(2010,2020,2030,2040)) +
    theme(legend.position='top')
}
plot.incidence = function(E){
  g = plot.out(E,'incidence',color='case.lab',scale=1000) +
    facet_wrap('.~groups',ncol=4,scales='free_y') +
    labs(x='Year',y='HIV Incidence (per 1000 PY)',color='',fill='') +
    theme(legend.position='top')
}
plot.obj.1 = function(X,y,y.lab,y.lims,t.lims=c(2005,2040),...){
  X$case.lab = factor(X$case.lab,labels=gsub('Left Behind: ','',levels(X$case.lab)))
  if (!is.null(t.lims)){ t = as.numeric(X$year); X = X[t>=t.lims[1] & t<=t.lims[2],] }
  g = ggplot(X[X$case.id!='bc',],aes_string(x='factor(year)',y=y,...)) +
    geom_boxplot(aes(color=case.lab),outlier.size=.5,lwd=.0,width=.6,position=position_dodge(.8)) +
    geom_boxplot(aes(fill=case.lab),outlier.color=NA,lwd=.3,width=.6,position=position_dodge(.8)) +
    labs(x='Year',y=y.lab,color='Left Behind:',fill='Left Behind:') + lims(y=y.lims) +
    scale_fill_manual(values=sget('cases.art','clr')) +
    scale_color_manual(values=sget('cases.art','clr')) +
    theme_light() +
    theme(legend.position=c(.01,.99),legend.justification=c(0,1))
  return(g)
}
# --------------------------------------------------------------------------------------------------
# numeric stuff - TODO: check
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
if (sys.nframe() == 0){
  E = load.expo.data()
  g = plot.incidence(E);                                    fig.save(uid,'ar_1_inc',w=12,h=4)
  g = plot.cascade(E,c('diagnosed','treated.c','vls.c'));   fig.save(uid,'art_1_cascade',w=8,h=6)
  g = plot.cascade(E,'vls.u',y.lab=spec$cascade$vls.u$lab); fig.save(uid,'art_1_vls',w=8,h=3)
  X = load.keyout.data()
  # numeric.obj.1(X) # TODO
  plot.obj.1(X,y='100*inf.red',y.lab='Infections averted (%)',      y.lims=c(0, 60)); fig.save(uid,'art_1_inf_red',w=8,h=4)
  plot.obj.1(X,y='100*inf.add',y.lab='Additional infections (%)',   y.lims=c(0,130)); fig.save(uid,'art_1_inf_add',w=8,h=4)
  plot.obj.1(X,y='100*inc.red',y.lab='Incidence reduction (%)',     y.lims=c(0,100)); fig.save(uid,'art_1_inc_red',w=8,h=4)
  plot.obj.1(X,y='  1*inc.add',y.lab='Additional incidence (times)',y.lims=c(0, 30)); fig.save(uid,'art_1_inc_add',w=8,h=4)
}