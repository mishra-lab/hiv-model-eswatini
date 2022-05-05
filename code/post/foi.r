source('post/config.r')

q.aes = function(q,g='case.lab',...){
  if (q == .5){
    return(aes_string(y='q0.5',color=g,linetype=g,...))
  } else {
    qq = function(op){ paste0(op,'(q',1-(1-q)/2,', q',(1-q)/2,')') }
    return(aes_string(ymin=qq('pmin'),ymax=qq('pmax'),fill=g,...))
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
  lab = list(foi='FOI Model',incidence='HIV Incidence (per person-year)',prevalence='HIV Prevalence')
  g = ggplot(X.,aes(x=t)) +
    geom_point(data=y.max,aes(x=t,y=q0.95),color=NA) +
    geom_ribbon(q.aes(.9),alpha=.25) +
    geom_line(q.aes(.5)) +
    scale_color_manual(values=sget('cases.foi','clr')) +
    scale_fill_manual(values=sget('cases.foi','clr')) +
    facet_wrap(.~facet,scales='free_y',drop=TRUE) +
    lims(x=c(1980,2030)) +
    labs(x='Year',y=lab[[out]],color=lab$foi,fill=lab$foi,linetype=lab$foi) +
    theme_light() + theme(legend.position='top')
  return(g)
}

X = read.csvs('foi-ep','expo','cases.foi')
pops = list(
  aq  = 'Non-Sex Work Women & Men',
  cli = 'Clients',
  fsw = 'Female Sex Workers')
X$facet = interaction(X$pop,X$op)
X$facet = factor(X$facet,
  levels=c('aq.raw','cli.raw','fsw.raw','aq.1-2','cli.1-2','fsw.1-2'),
  labels=c(pops$aq,pops$cli,pops$fsw,paste0(c('',' ','  '),'Difference vs <4*>',c('',' ','  '))))
plot.ep(X,'incidence',pops); fig.save(uid,paste0('foi.ep.incidence'),w=10,h=7)
plot.ep(X,'prevalence',pops); fig.save(uid,paste0('foi.ep.prevalence'),w=10,h=7)
