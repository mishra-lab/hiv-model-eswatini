library('reshape2')
source('post/config.r')
source('post/wiw.r')

out.labs = list(
  cuminfect = 'Cumulative Additional Infections',
  incidence = 'Additional Incidence',
  diagnosed = 'Diagnosed among PLHIV',
  treated_c = 'Treated among Diagnosed',
  vls_c     = 'VLS among Treated',
  treated_u = 'Treated among PLHIV',
  vls_u     = 'VLS among PLHIV')

pop.labs = list(
  all = 'Overall',
  w   = 'Women Overall',
  m   = 'Men Overall',
  aq  = 'Lower Risk',
  fsw = 'FSW',
  cli = 'Clients')

melt.expo.cascade = function(X.wide,...){
  X = melt.expo.s(X.wide,pop=c('all','aq','fsw','cli'),
    out=c('diagnosed','treated_c','vls_c','treated_u','vls_u'),...)
  X$value = 100 * X$value
  X$out = factor(X$out,names(out.labs),out.labs)
  X$pop = factor(X$pop,names(pop.labs),pop.labs)
  return(X)
}

plot.clean.art = function(g){
  g = g + labs(x='Year',color='Left behind:',fill='Left behind:',linetype='Left behind:') +
    scale_color_manual(values=set.cols$art) +
    scale_fill_manual(values=set.cols$art)
  return(g)
}

plot.1.rai = function(X){
  X$out = factor(X$out,names(out.labs),out.labs)
  g = plot.expo.box(expo.qs(X),unlist(out.labs),'case.lab',seq(2005,2030,5)) +
    facet_grid('out',scales='free_y') + 
    labs(y='Relative Additional Infections (%)') +
    theme(legend.pos=c(.005,.99),legend.justification=c(0,1))
  g = plot.clean.art(g)
  fig.save(uid,N$sam,'art.1.rai',w=6,h=5)
}

num.1.rai = function(X,t.hor=2030){
  X.hor = X[X$t==t.hor,]
  print(expo.qs(X.hor,q=c(.025,.5,.975)))
  x.case = function(case){ X.hor[X.hor$case==case,] }
  X.ref = x.case('fsw-cli-')
  X.ref$value = (X.ref$value - x.case('fsw+cli+')$value) / X.ref$value
  print(expo.qs(X.ref,q=c(.025,.5,.975)))
}

main.1.rai = function(){
  X.wide = read.csvs('art-rf','expo','art',rdata='load')
  X = melt.expo.s(X.wide,pop='all',out=c('incidence','cuminfect'))
  base.value = X[X$case=='base',]$value
  X$value = 100 * (X$value - base.value) / base.value
  X = X[X$case!='base',]
  num.1.rai(X)
  plot.1.rai(X)
}

main.1.wiw = function(){
  X = clean.wiw.data(read.csvs('art-rf','wiw-diff','art',skip='base'))
  X = X[X$t > 2003 & X$t <= 2030,]
  for (margin in c('part','from','to')){
    g = do.margin(X,margin,type='abs',strat='case.lab') +
      facet_grid('~case.lab',scales='free_y') +
      labs(y='Additional Infections (\'000s)') +
      guides(colour=guide_legend(nrow=1),fill=guide_legend(nrow=1))
    fig.save(uid,N$sam,'art.wiw',margin, w=8,h=3)
  }
}

main.1.cascade = function(){
  X.wide = read.csvs('art-rf','expo','art',rdata='load')
  X.wide = X.wide[X.wide$t > 1995 & X.wide$t <= 2030,]
  X = melt.expo.cascade(X.wide)
  g = plot.expo.ribbon(expo.qs(X),unlist(out.labs),'case.lab') +
    facet_grid('out ~ pop') + lims(y=c(0,100)) +
    labs(y='Cascade Step (%)') +
    theme(legend.position='top')
  g = plot.clean.art(g)
  fig.save(uid,N$sam,'art.1.cascade',w=8,h=10)
}

main.2.cascade = function(){
  X.wide = read.csvs('art-ss','expo','sens',rdata='load')
  X = melt.expo.cascade(X.wide,t=2020)
  g = ggplot(X,aes(x=value,y=pop,color=pop,fill=pop)) +
    geom_boxplot(alpha=.25,outlier.alpha=1,outlier.shape=3,outlier.size=.5) +
    facet_grid('out') +
    scale_x_continuous(breaks=seq(0,100,20),limits=c(0,100)) +
    scale_color_manual(values=set.cols$pop.art) +
    scale_fill_manual(values=set.cols$pop.art) +
    labs(x='Cascade step (%)',y='Population',color='',fill='')
  g = plot.clean(g,legend.position='none')
  fig.save(uid,N$sam,'art.2.cascade',w=5,h=9)
}

# main.1.cascade()
# main.1.rai()
# main.1.wiw()
# main.2.cascade()
