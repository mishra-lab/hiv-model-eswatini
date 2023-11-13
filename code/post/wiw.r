clean.wiw.data = function(X,infs='q0.5'){
  X$infs = X[[infs]]
  X = X[,!grepl('^q\\d',names(X))]
  X$from  = factor(1+X$fs*4+X$fi,levels=1:8,labels=set.labs$pop.all)
  X$to    = factor(1+X$ts*4+X$ti,levels=1:8,labels=set.labs$pop.all)
  X$part  = factor(1+X$p,levels=1:4,labels=set.labs$part)
  return(X)
}
do.norm = function(X,margin,strat=NULL,scale=100){
  X.s = aggregate(formula(paste(c('infs ~ t',strat),collapse=' + ')),X,sum)
  X = merge(X,rename.cols(X.s,infs='infs.total'))
  X$infs = scale * X$infs / X$infs.total
  X$infs.total = NULL
  return(X)
}
do.alluvial = function(X,fill='part',norm=FALSE){
  # paste hack to maintain order in alluvial
  X$from = factor(1+X$fs*4+X$fi,1:8,paste(set.labs$pop.all,''))
  X$to   = factor(1+(1-X$ts)*4+X$ti,1:8,paste0('',c(set.labs$pop.all[5:8],set.labs$pop.all[1:4])))
  X.long = to_lodes_form(X,key='pop',axis=c('from','to'))
  if (norm){ X.long = do.norm(X.long,margin=fill,scale=200) }
  ylab = paste('Yearly Infections',ifelse(norm,'(%)','(\'000s)'))
  g = ggplot(X.long,aes(y=infs,x=pop,stratum=stratum,label=stratum,alluvium=alluvium)) +
    geom_alluvium(aes_string(fill=fill),width=.4) +
    geom_stratum(color=rgb(1,1,1),alpha=.5,width=.4) +
    geom_text(stat='stratum',size=2.5) +
    labs(x='Population',y=ylab,fill='Partnership Type') +
    scale_x_discrete(expand=c(.12,.12)) +
    scale_y_continuous(expand=c(.02,.02)) +
    scale_fill_manual(values=set.cols[[fill]])
  g = plot.clean(g)
}
do.alluvial.facet = function(X,tvec=seq(1990,2020,10),nrow=1){
  g = do.alluvial(X[X$t %in% tvec,],norm=TRUE) +
    facet_wrap(~factor(t),nrow=nrow) +
    theme(legend.position='top')
}
do.ratio = function(X){
  X.p = aggregate(infs~t+from+to,X,sum) # sum parts
  X. = merge(rename.cols(aggregate(infs~t+from,X.p,sum),infs='infs.fr',from='pop'),
             rename.cols(aggregate(infs~t+to,  X.p,sum),infs='infs.to',to='pop'))
  X.$infs.ratio = X.$infs.fr / X.$infs.to
  b = 2^seq(-6,2)
  g = ggplot(X.,aes(y=infs.ratio,x=t,color=pop)) +
    geom_hline(yintercept=1,color=rgb(.8,.8,.8),lwd=1) +
    geom_line(lwd=.7) +
    scale_y_continuous(trans='log2',breaks=b,labels=MASS::fractions(b),lim=range(b)) +
    scale_color_manual(values=set.cols$pop.all) +
    labs(x='Year',y='Yearly Infections Transmitted / Acquired ',color='',fill='')
  g = plot.clean(g,legend.position='top')
}
do.margin = function(X,margin,type='both',strat=NULL){
  X.abs = aggregate(formula(paste(c('infs ~ t',margin,strat),collapse=' + ')),X,sum)
  X.rel = do.norm(X.abs,margin=margin,strat=strat)
  X.abs$scale = 'Absolute (\'000s)'
  X.rel$scale = 'Proportion (%)'
  if (margin %in% c('from','to')){
    i.wom = grepl('Women\\s|\\sFSW',X.abs[[margin]])
    X.abs$scale[i.wom]  = paste0(X.abs$scale[i.wom],', Women')
    X.abs$scale[!i.wom] = paste0(X.abs$scale[!i.wom],', Men')
    clr = 'pop.all'
  } else { clr = margin }
  X.long = switch(type,abs=X.abs,rel=X.rel,both=rbind(X.abs,X.rel))
  g = ggplot(X.long,aes_string(x='t',y='infs',color=margin,fill=margin)) +
    geom_area(alpha=.7,size=.1) +
    scale_color_manual(values=set.cols[[clr]]) +
    scale_fill_manual(values=set.cols[[clr]]) +
    facet_grid('scale~.',scales='free_y') +
    labs(x='Year',y='Yearly Infections',color='',fill='')
  g = plot.clean(g,legend.position='top')
}
