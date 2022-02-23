library('gifski')
source('post/config.r')

# --------------------------------------------------------------------------------------------------
load.data = function(diff=FALSE){
  if (diff){
    X = load.csvs('fit','infs-diff',cases[1:4])
  } else {
    X = load.csvs('fit','infs')
  }
  # X = load.csvs('sens','infs',cases.sens[1],b=seq(N$batch)-1) # DEBUG: obj 2
  # X = aggregate(infections~t+fs+fi+ts+ti+p,X,median) # DEBUG: obj 2
  X = X[X$infections>0,] # for speed
  X$from  = factor(1+X$fs*4+X$fi,levels=1:8,labels=lab$pop)
  X$to    = factor(1+X$ts*4+X$ti,levels=1:8,labels=lab$pop)
  X$part  = factor(1+X$p,levels=1:4,labels=lab$part)
  return(X)
}
do.norm = function(X,n,scale=100){
  total.inf = aggregate(infections~t,X,sum)$infections
  X$infections = scale * X$infections / rep(total.inf,n)
  return(X)
}
alluvial = function(X,fill='part',norm=FALSE){
  # paste hack to maintain order in alluvial
  X$from = factor(1+X$fs*4+X$fi,levels=1:8,labels=paste(lab$pop,''))
  X$to   = factor(1+(1-X$ts)*4+X$ti,levels=1:8,labels=paste0('',c(lab$pop[5:8],lab$pop[1:4])))
  X.long = to_lodes_form(X,key='pop',axis=c('from','to'))
  if (norm){ X.long = do.norm(X.long,sum(X$t==X$t[1]),200) }
  ylab = paste('Yearly Infections',ifelse(norm,'(%)','(\'000s)'))
  g = ggplot(X.long,aes(y=infections,x=pop,stratum=stratum,label=stratum,alluvium=alluvium)) +
    geom_alluvium(aes_string(fill=fill),width=.4) +
    geom_stratum(color=rgb(.5,.5,.5),alpha=.5,width=.4) +
    geom_text(stat='stratum',size=3) +
    labs(x='Population',y=ylab,fill='') +
    scale_x_discrete(expand=c(.15,.15)) +
    scale_fill_manual(values=clr[[fill]]) +
    theme_light()
  return(g)
}
do.alluvial.facet = function(X){
  tvec = seq(1990,2040,10)
  g = alluvial(X[X$t %in% tvec,],norm=TRUE) +
    facet_wrap(~factor(t)) +
    theme(legend.position='top')
}
do.alluvial.gif = function(X){
  save_gif(lapply(unique(X$t),function(t){ # print(t)
    print(alluvial(X[X$t==t,],norm=TRUE) + ggtitle(t))
  }),'alluvial.gif',width=600,height=800,delay=1/10,loop=TRUE)
}
do.ratio = function(X){
  X.p = aggregate(infections~t+from+to,X,sum) # sum parts
  X. = merge(rename.cols(aggregate(infections~t+from,X.p,sum),infections='inf.fr',from='pop'),
             rename.cols(aggregate(infections~t+to,  X.p,sum),infections='inf.to',to='pop'))
  X.$inf.ratio = X.$inf.fr / X.$inf.to
  g = ggplot(X.,aes(y=inf.ratio,x=t,color=pop)) +
    geom_hline(yintercept=1,color=rgb(.8,.8,.8),lwd=1) +
    geom_line(lwd=.7) +
    scale_color_manual(values=clr[['from']]) +
    labs(x='Year',y='Yearly Infections Transmitted / Acquired (Ratio)',color='',fill='') +
    theme_light()
  return(g)
}
do.margin = function(X,margin,rel=TRUE){
  # partnership types
  X.abs = aggregate(formula(paste('infections ~ t + case.lab +',margin)),X,sum)
  X.rel = do.norm(X.abs,sum(X.abs$t==X.abs$t[1]))
  X.abs$scale = 'Absolute (\'000s)'
  X.rel$scale = 'Proportion (%)'
  if (margin == 'part'){ print(X.rel[X.rel$t %in% seq(1990,2040,5),]) }
  if (margin %in% c('from','to')){
    i.men = grepl('Men\\s|\\sClients',X.abs[[margin]])
    X.abs$infections[i.men] = -X.abs$infections[i.men]
  }
  if (rel){ X.long = rbind(X.abs,X.rel) } else { X.long = X.abs }
  g = ggplot(X.long,aes_string(x='t',y='infections',color=margin,fill=margin)) +
    geom_area(alpha=.7,size=.1) +
    scale_color_manual(values=clr[[margin]]) +
    scale_fill_manual(values=clr[[margin]]) +
    labs(x='Year',y='Infections',color='',fill='') +
    theme_light()
  if (margin %in% c('from','to')){ g = g + scale_y_continuous(labels=abs) }
  if (rel){ g = g + facet_wrap(vars(scale),scales='free_y') }
  return(g)
}

# --------------------------------------------------------------------------------------------------
X = load.data()
X.base = X[X$case=='base',]
do.ratio(X.base);          fig.save(uid,'inf-ratio',w=5,h=4)
do.margin(X.base,'part');  fig.save(uid,'inf-part',w=8,h=4)
do.margin(X.base,'to');    fig.save(uid,'inf-to'  ,w=8,h=4)
do.margin(X.base,'from');  fig.save(uid,'inf-from',w=8,h=4)
do.alluvial.facet(X.base); fig.save(uid,'inf-alluvial',w=12,h=16)
# do.alluvial.gif(X.base)
X.diff = load.data(diff=TRUE)
do.margin(X.diff,'part',rel=FALSE) + facet_grid(cols=vars(case.lab)); fig.save(uid,'inf-diff-part',w=10,h=3)
do.margin(X.diff,'to',  rel=FALSE) + facet_grid(cols=vars(case.lab)); fig.save(uid,'inf-diff-to',  w=10,h=3)
do.margin(X.diff,'from',rel=FALSE) + facet_grid(cols=vars(case.lab)); fig.save(uid,'inf-diff-from',w=10,h=3)
