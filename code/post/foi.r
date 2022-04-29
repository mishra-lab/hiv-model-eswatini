source('post/config.r')

q.aes = function(q,g='case.lab',...){
  if (q == .5){
    return(aes_string(y='q0.5',color=g,linetype=g,...))
  } else {
    qq = function(op){ paste0(op,'(q',1-(1-q)/2,', q',(1-q)/2,')') }
    return(aes_string(ymin=qq('pmin'),ymax=qq('pmax'),fill=g,...))
  }
}
plot.ep = function(X,out,pops){
  X. = X[X$out==out,]
  X.$pop.lab = factor(X.$pop,levels=names(pops),labels=unname(pops))
  X. = X.[!is.na(X.$pop.lab),]
  lab = list(foi='FOI Model',incidence='HIV Incidence (per person-year)',prevalence='HIV Prevalence')
  g = ggplot(X.,aes(x=t)) +
    geom_ribbon(q.aes(.9),alpha=.25) +
    geom_line(q.aes(.5)) +
    scale_color_manual(values=sget('cases.foi','clr')) +
    scale_fill_manual(values=sget('cases.foi','clr')) +
    facet_wrap(.~facet,scales='free_y',drop=TRUE) +
    lims(x=c(1980,2030)) +
    labs(x='Year',y=lab[[out]],color=lab$foi,fill=lab$foi,linetype=lab$foi) +
    theme_light()
  return(g)
}

X = load.csvs('foi','ep','cases.foi')
pops = list(
  aq  = 'Non-Sex Work Women & Men',
  cli = 'Clients',
  fsw = 'Female Sex Workers')
X$facet = interaction(X$pop,X$op)
X$facet = factor(X$facet,
  levels=c('aq.raw','cli.raw','fsw.raw','aq.1-2','cli.1-2','fsw.1-2'),
  labels=c(pops$aq,pops$cli,pops$fsw,'Difference',' Difference ','  Difference  '))
plot.ep(X,'incidence',pops); fig.save(uid,paste0('ep-incidence'),w=12,h=7)
