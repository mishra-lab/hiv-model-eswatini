source('post/config.r')

lab = list(
  incidence  = 'HIV Incidence (per person-year)',
  prevalence = 'HIV Prevalence',
  msp   = 'Main / Spousal',
  cas   = 'Casual',
  swq   = 'All Sex Work',
  aqfr  = 'Non-Sex Work Women & Men', aqto  = 'Non-Sex Work Women & Men',
  fswfr = 'Female Sex Workers',       fswto = 'Female Sex Workers',
  clifr = 'Clients',                  clito = 'Clients',
  foi = 'Approach (median, 90% confidence interval)')

q.aes = function(q,grp='case.lab',t='t',...){
  if (q == 0){
    return(aes_string(y='q0.5',color=grp,linetype=grp,...))
  } else {
    qq = function(qi,op){ paste0(op,'(q',1-(1-qi)/2,', q',(1-qi)/2,')') }
    if (is.numeric(q)) {
      return(aes_string(ymin=qq(q,'pmin'),ymax=qq(q,'pmax'),x=t,fill=grp,...))
    } else {
      return(aes_string(ymin=qq(.95,'pmin'),ymax=qq(.95,'pmax'),fill=grp,color=grp,
                        lower=qq(.5,'pmin'),upper=qq(.5,'pmax'),middle='q0.5',
                        x=paste0('factor(',t,')'),group=paste0('interaction(',t,',',grp,')'),...))
    }
  }
}
make.y.max = function(X.,q,m=.55){
  ystr = paste0('q',1-(1-q)/2)
  y.max = aggregate(formula(paste(ystr,'~pop')),X.[X.$op=='raw',],max)
  y.max$t = 2000
  y.max$facet = unique(X.$facet[X.$op=='1-2/2'])
  y.max = rbind(y.max,y.max)
  y.max[[ystr]] = y.max[[ystr]] * rep(c(-m,m),each=nrow(y.max)/2)
  return(y.max)
}
plot.expo = function(X,out,pops,grp='case.lab',clrs,lts,diff=TRUE){
  if (missing(clrs)){ clrs = sget('cases.foi','clr') }
  if (missing(lts)){ lts = sget('cases.foi','lt') }
  X. = X[X$out==out,]
  X.$pop.lab = factor(X.$pop,levels=names(pops),labels=unname(pops))
  X. = X.[!is.na(X.$pop.lab),]
  g = ggplot(X.,aes(x=t)) +
    geom_ribbon(q.aes(.9,grp=grp),alpha=.25) +
    geom_line(q.aes(0,grp=grp)) +
    scale_color_manual(values=clrs) +
    scale_fill_manual(values=clrs) +
    scale_linetype_manual(values=lts) +
    facet_wrap(.~facet,scales='free_y',drop=TRUE) +
    lims(x=c(1980,2030)) +
    labs(x='Year',y=lab[[out]],color=lab$foi,fill=lab$foi,linetype=lab$foi) +
    theme_light() + theme(legend.position='top')
  if (diff){ g = g + geom_point(data=make.y.max(X.,.9),aes(x=t,y=q0.95),color=NA) }
  return(g)
}
plot.tpaf = function(X,t,tpaf.pop,tpaf.t0,grp='case.lab',inc.pop='all',mode='box',lp='top',t.type='t.hor',clrs){
  if (missing(t)) { t = unique(X[[t.type]]) }
  if (missing(clrs)){ clrs = sget('cases.foi','clr') }
  tlab = list(t='Year',t.hor='Time Horizon (Years)')
  X. = X[(X$tpaf.pop %in% tpaf.pop) & (X$tpaf.t0 %in% tpaf.t0) & (X$pop == inc.pop) & (X[[t.type]] %in% t),]
  X.$tpaf.pop = factor(X.$tpaf.pop,levels=tpaf.pop,labels=lab[tpaf.pop])
  g = ggplot(X.) +
    scale_color_manual(values=clrs) +
    scale_fill_manual(values=clrs) +
    labs(x=tlab[[t.type]],y='tPAF',fill=lab$foi,color=lab$foi,linetype=lab$foi) +
    theme_light() + theme(legend.position=lp)
  if (length(tpaf.pop)  > 1 & length(tpaf.t0)  > 1){ g = g + facet_grid(tpaf.t0~tpaf.pop) }
  if (length(tpaf.pop)  > 1 & length(tpaf.t0) == 1){ g = g + facet_grid(.~tpaf.pop) }
  if (length(tpaf.pop) == 1 & length(tpaf.t0)  > 1){ g = g + facet_grid(.~tpaf.t0) }
  if (mode=='box'){ g = g + geom_boxplot(q.aes('box',grp=grp,t=t.type),stat='identity',width=.7,alpha=.25,pos=position_dodge2(pad=.3)) }
  if (mode=='rib'){ g = g + geom_ribbon(q.aes(.9,grp=grp,t=t.type),alpha=.25) + geom_line(q.aes(0,grp=grp,x=t.type)) }
  return(g)
}
main.expo = function(param='ep'){
  # param: ['ep','fit']
  diff = switch(param,ep=TRUE,fit=FALSE)
  X = read.csvs(paste0('foi-',param),'expo','cases.foi')
  pops = list(
    aq  = 'Non-Sex Work Women & Men',
    cli = 'Clients',
    fsw = 'Female Sex Workers')
  X$facet = interaction(X$pop,X$op)
  X$facet = factor(X$facet,
    levels=c('aq.raw','cli.raw','fsw.raw','aq.1-2','cli.1-2','fsw.1-2'),
    labels=c(pops$aq,pops$cli,pops$fsw,paste0(c('',' ','  '),'Difference vs <4*>',c('',' ','  '))))
  plot.expo(X,'incidence',pops,diff=diff); fig.save(uid,paste0('foi.',param,'.incidence'),w=10,h=4+3*diff)
  plot.expo(X,'prevalence',pops,diff=diff); fig.save(uid,paste0('foi.',param,'.prevalence'),w=10,h=4+3*diff)
}
main.smdm = function(){
  ep.smdm = function(){
    X = read.csvs('foi-ep','expo','cases.foi',skip=c('bmy','bpd'))
    X = X[X$op!='1-2',]
    X$model = factor(interaction(X$case.id,X$op),levels=c('bpy.raw','base.raw','bpy.1-2/2'),
      labels=c('Instantaneous','Proposed','Relative Difference  '))
    X$facet = factor(interaction(X$pop,X$op),levels=c('all.raw','all.1-2/2'),
      labels=c('Absolute','Instantaneous - Proposed / Proposed'))
    clrs = c(rgb(.267,.005,.329),rgb(.993,.906,.144),rgb(.736,.216,.330))
    plot.expo(X,'incidence',list(all='Overall'),g='model',clrs=clrs) +
      scale_linetype_manual(values=c('32','solid','62'));
    fig.save(uid,'foi.ep.incidence.smdm',w=7,h=4)
  }
  tpaf.smdm = function(){
    X = read.csvs('foi-tpaf','expo','cases.foi',skip=c('bmy','bpd'))
    X$model = factor(X$case.id,levels=c('base','bpy'),labels=c('Proposed','Instantaneous'))
    X$t.hor = X$t - X$tpaf.t0
    t = c(1,2,5,10)
    clrs = c(rgb(.978,.558,.035),rgb(.342,.062,.429))
    g = plot.tpaf(X,t=t,tpaf.pop=c('msp','cas','swq'),tpaf.t0=2000,g='model',clrs=clrs);
    fig.save(uid,'foi.tpaf.parts.smdm',w=6,h=3.5)
  }
  ep.smdm()
  # tpaf.smdm()
}
main.tpaf = function(){
  t0 = c(1990,2000,2010)
  t = c(1,2,5,10)
  X = read.csvs('foi-tpaf','expo','cases.foi')
  X$t.hor = X$t - X$tpaf.t0
  g = plot.tpaf(X,t=t,tpaf.pop=c('msp','cas','swq'),     tpaf.t0=2000);           fig.save(uid,'foi.tpaf.3.parts',w=10,h=3.5)
  g = plot.tpaf(X,t=t,tpaf.pop=c('aqfr','clifr','fswfr'),tpaf.t0=2000,lp='none'); fig.save(uid,'foi.tpaf.3.popfr',w=10,h=3)
  g = plot.tpaf(X,t=t,tpaf.pop=c('aqto','clito','fswto'),tpaf.t0=2000,lp='none'); fig.save(uid,'foi.tpaf.3.popto',w=10,h=3)
  g = plot.tpaf(X,t=t,tpaf.pop=c('aqfr','clifr','fswfr'),tpaf.t0=t0); fig.save(uid,'foi.tpaf.33.popfr',w=10,h=8)
  g = plot.tpaf(X,t=t,tpaf.pop=c('aqto','clito','fswto'),tpaf.t0=t0); fig.save(uid,'foi.tpaf.33.popto',w=10,h=8)
  g = plot.tpaf(X,t=t,tpaf.pop=c('msp','cas','swq'),     tpaf.t0=t0); fig.save(uid,'foi.tpaf.33.parts',w=10,h=8)
}
main.smdm()
main.expo('ep')
main.expo('fit')
main.tpaf()