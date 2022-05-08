source('post/config.r')

lab = list(
  incidence  = 'HIV Incidence (per person-year)',
  prevalence = 'HIV Prevalence',
  msp   = 'Main / Spousal',
  cas   = 'Casual',
  swq   = 'Sex Work',
  aqfr  = 'Non-Sex Work Women & Men', aqto  = 'Non-Sex Work Women & Men',
  fswfr = 'Female Sex Workers',       fswto = 'Female Sex Workers',
  clifr = 'Clients',                  clito = 'Clients',
  foi = 'FOI Model')

q.aes = function(q,g='case.lab',t='t',...){
  if (q == 0){
    return(aes_string(y='q0.5',color=g,linetype=g,...))
  } else {
    qq = function(qi,op){ paste0(op,'(q',1-(1-qi)/2,', q',(1-qi)/2,')') }
    if (is.numeric(q)) {
      return(aes_string(ymin=qq(q,'pmin'),ymax=qq(q,'pmax'),x=t,fill=g,...))
    } else {
      return(aes_string(ymin=qq(.95,'pmin'),ymax=qq(.95,'pmax'),fill=g,color=g,
                        lower=qq(.5,'pmin'),upper=qq(.5,'pmax'),middle='q0.5',
                        x=paste0('factor(',t,')'),group=paste0('interaction(',t,',',g,')'),...))
    }
  }
}
make.y.max = function(X.,q,m=.55){
  ystr = paste0('q',1-(1-q)/2)
  y.max = aggregate(formula(paste(ystr,'~pop')),X.[X.$op=='raw',],max)
  y.max$t = 2000
  y.max$facet = unique(X.$facet[X.$op=='1-2'])
  y.max = rbind(y.max,y.max)
  y.max[[ystr]] = y.max[[ystr]] * rep(c(-m,m),each=nrow(y.max)/2)
  return(y.max)
}
plot.ep = function(X,out,pops){
  X. = X[X$out==out,]
  X.$pop.lab = factor(X.$pop,levels=names(pops),labels=unname(pops))
  X. = X.[!is.na(X.$pop.lab),]
  y.max = make.y.max(X.,.9)
  g = ggplot(X.,aes(x=t)) +
    geom_point(data=y.max,aes(x=t,y=q0.95),color=NA) +
    geom_ribbon(q.aes(.9),alpha=.25) +
    geom_line(q.aes(0)) +
    scale_color_manual(values=sget('cases.foi','clr')) +
    scale_fill_manual(values=sget('cases.foi','clr')) +
    facet_wrap(.~facet,scales='free_y',drop=TRUE) +
    lims(x=c(1980,2030)) +
    labs(x='Year',y=lab[[out]],color=lab$foi,fill=lab$foi,linetype=lab$foi) +
    theme_light() + theme(legend.position='top')
  return(g)
}
main.ep = function(){
  X = read.csvs('foi-ep','expo','cases.foi')
  pops = list(
    aq  = 'Non-Sex Work Women & Men',
    cli = 'Clients',
    fsw = 'Female Sex Workers')
  X$facet = interaction(X$pop,X$op)
  X$facet = factor(X$facet,
    levels=c('aq.raw','cli.raw','fsw.raw','aq.1-2','cli.1-2','fsw.1-2'),
    labels=c(pops$aq,pops$cli,pops$fsw,paste0(c('',' ','  '),'Difference vs <4*>',c('',' ','  '))))
  plot.ep(X,'incidence',pops); fig.save(uid,'foi.ep.incidence',w=10,h=7)
  plot.ep(X,'prevalence',pops); fig.save(uid,'foi.ep.prevalence',w=10,h=7)
}
plot.tpaf = function(X,t,tpaf.pop,tpaf.t0,inc.pop='all',mode='box',lp='top',t.type='t.hor'){
  if (missing(t)) { t = unique(X[[t.type]]) }
  tlab = list(t='Year',t.hor='Time Horizon (Years)')
  X. = X[(X$tpaf.pop %in% tpaf.pop) & (X$tpaf.t0 %in% tpaf.t0) & (X$pop == inc.pop) & (X[[t.type]] %in% t),]
  X.$tpaf.pop = factor(X.$tpaf.pop,levels=tpaf.pop,labels=lab[tpaf.pop])
  g = ggplot(X.) +
    scale_color_manual(values=sget('cases.foi','clr')) +
    scale_fill_manual(values=sget('cases.foi','clr')) +
    labs(x=tlab[[t.type]],y='tPAF',fill=lab$foi,color=lab$foi,linetype=lab$foi) +
    theme_light() + theme(legend.position=lp)
  if (length(tpaf.pop)  > 1 & length(tpaf.t0)  > 1){ g = g + facet_grid(tpaf.t0~tpaf.pop) }
  if (length(tpaf.pop)  > 1 & length(tpaf.t0) == 1){ g = g + facet_grid(.~tpaf.pop) }
  if (length(tpaf.pop) == 1 & length(tpaf.t0)  > 1){ g = g + facet_grid(.~tpaf.t0) }
  if (mode=='box'){ g = g + geom_boxplot(q.aes('box',t=t.type),stat='identity',width=.7,alpha=.25,pos=position_dodge2(pad=.3)) }
  if (mode=='rib'){ g = g + geom_ribbon(q.aes(.9,t=t.type),alpha=.25) + geom_line(q.aes(0,x=t.type)) }
  return(g)
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
main.ep()
main.tpaf()