library('gifski')
source('post/config.r')

lab = list(pop=sget('pops','lab'),part=sget('parts','lab'))
clr = list(from=sget('pops','clr'),to=sget('pops','clr'),part=sget('parts','clr'))
# --------------------------------------------------------------------------------------------------
clean.data = function(X){
  # X = load.csvs('sens','infs',cases.sens[1],b=seq(N$batch)-1) # DEBUG: obj 2
  # X = aggregate(infections~t+fs+fi+ts+ti+p,X,median) # DEBUG: obj 2
  X = X[X$infections>0,] # for speed
  X$from  = factor(1+X$fs*4+X$fi,levels=1:8,labels=lab$pop)
  X$to    = factor(1+X$ts*4+X$ti,levels=1:8,labels=lab$pop)
  X$part  = factor(1+X$p,levels=1:4,labels=lab$part)
  return(X)
}
do.norm = function(X,margin,strat='1',scale=100){
  if (strat != '1'){
    X = X[order(X[[margin]],X[[strat]],X$t),]
  }
  total.inf = aggregate(formula(paste('infections ~ t +',strat)),X,sum)$infections
  X$infections = scale * X$infections / total.inf
  return(X)
}
alluvial = function(X,fill='part',norm=FALSE){
  # paste hack to maintain order in alluvial
  X$from = factor(1+X$fs*4+X$fi,levels=1:8,labels=paste(lab$pop,''))
  X$to   = factor(1+(1-X$ts)*4+X$ti,levels=1:8,labels=paste0('',c(lab$pop[5:8],lab$pop[1:4])))
  X.long = to_lodes_form(X,key='pop',axis=c('from','to'))
  if (norm){ X.long = do.norm(X.long,margin=fill,scale=200) }
  ylab = paste('Yearly Infections',ifelse(norm,'(%)','(\'000s)'))
  g = ggplot(X.long,aes(y=infections,x=pop,stratum=stratum,label=stratum,alluvium=alluvium)) +
    geom_alluvium(aes_string(fill=fill),width=.4) +
    geom_stratum(color=rgb(1,1,1),alpha=.5,width=.4) +
    geom_text(stat='stratum',size=2.5) +
    labs(x='Population',y=ylab,fill='Partnership Type') +
    scale_x_discrete(expand=c(.12,.12)) +
    scale_y_continuous(expand=c(.02,.02)) +
    scale_fill_manual(values=clr[[fill]]) +
    theme_light()
  return(g)
}
do.alluvial.facet = function(X,tvec=seq(1990,2040,10)){
  g = alluvial(X[X$t %in% tvec,],norm=TRUE) +
    facet_wrap(~factor(t)) +
    theme(legend.position='bottom',
      strip.text=element_text(color='black'),
      strip.background=element_rect(fill='white'))
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
    scale_y_continuous(trans='log10',breaks=c(.1,.2,.5,1,2,5,10)) +
    labs(x='Year',y='Yearly Infections Transmitted / Acquired',color='',fill='') +
    theme_light()
  return(g)
}
do.margin = function(X,margin,type='both',strat='case.lab'){
  # partnership types
  X.abs = aggregate(formula(paste('infections ~ t + ',margin,' + ',strat)),X,sum)
  X.rel = do.norm(X.abs,margin=margin,strat=strat)
  X.abs$scale = 'Absolute (\'000s)'
  X.rel$scale = 'Proportion (%)'
  if (margin == 'part'){ print(X.rel[X.rel$t %in% seq(1990,2040,5),]) }
  if (margin == 'part'){ print(aggregate(infections~part+case.lab,X.rel,max)) }
  if (margin %in% c('from','to')){
    i.men = grepl('Men\\s|\\sClients',X.abs[[margin]])
    X.abs$infections[i.men] = -X.abs$infections[i.men]
  }
  X.long = switch(type,abs=X.abs,rel=X.rel,both=rbind(X.abs,X.rel))
  ylab   = switch(type,abs='(\'000s)',rel='(%)','')
  g = ggplot(X.long,aes_string(x='t',y='infections',color=margin,fill=margin)) +
    geom_area(alpha=.7,size=.1) +
    scale_color_manual(values=clr[[margin]]) +
    scale_fill_manual(values=clr[[margin]]) +
    labs(x='Year',y=paste('Infections',ylab),color='',fill='') +
    theme_light()
  if (margin %in% c('from','to')){ g = g + scale_y_continuous(labels=abs) }
  if (type=='both'){ g = g + facet_wrap(vars(scale),scales='free_y') }
  return(g)
}
# --------------------------------------------------------------------------------------------------
main.ims.poster = function(){
  X = clean.data(read.csv(gen.name('fit','infs','base','all')))
  X.base = X[X$t<=2025,]
  do.alluvial.facet(X.base,c(1990,2005,2020)); fig.save(uid,'inf-base-alluvial-3',w=9,h=6)
  do.margin(X.base,'part','rel') + guides(color='none',fill='none'); fig.save(uid,'inf-base-part-rel',w=3,h=3)
  do.margin(X.base,'from','rel') + guides(color='none',fill='none'); fig.save(uid,'inf-base-from-rel',w=3,h=3)
  do.margin(X.base,'to','rel');                                      fig.save(uid,'inf-base-to-rel',w=4.3,h=3)
}
main.foi = function(){
  X = clean.data(read.csvs('foi-fit','infs','cases.foi'))
  X = X[X$t>=1985 & X$t<=2030,]
  X.base = X[X$case=='base',]
  g.clean = function(g,lp){ g + facet_grid(cols=vars(case.lab)) + theme(legend.position=lp) }
  g.clean(do.margin(X,'part','rel','case.lab'),'top');  fig.save(uid,'inf-foi-part-rel',w=10,h=4)
  g.clean(do.margin(X,'from','rel','case.lab'),'top');  fig.save(uid,'inf-foi-from-rel',w=10,h=4)
  g.clean(do.margin(X,'to',  'rel','case.lab'),'none'); fig.save(uid,'inf-foi-to-rel',w=10,h=3.2)
  do.alluvial.facet(X.base,c(1990,2005,2020)) + theme(legend.position='top'); fig.save(uid,'inf-base-alluvial-3',w=9,h=7)
  # TODO: base without rel?
}
main.art = function(){
  # infs
  X = clean.data(read.csvs('art','infs','cases.art'))
  X = X[X$t>=1985 & X$t<=2040,]
  X.base = X[X$case=='base',]
  do.margin(X.base,'part'); fig.save(uid,'inf-base-part-both',w=7,h=3);
  q()
  do.margin(X.base,'from'); fig.save(uid,'inf-base-from-both',w=7,h=3);
  do.margin(X.base,'to');   fig.save(uid,'inf-base-to-both',w=7,h=3);
  do.ratio(X.base);         fig.save(uid,'inf-base-ratio',w=5,h=3.5)
  do.alluvial.facet(X.base,seq(1990,2040,10)) + theme(legend.position='top'); fig.save(uid,'inf-base-alluvial',w=8,h=10)
  # infs-diff
  X.diff = clean.data(read.csvs('art','infs-diff','cases.art',skip='base'))
  X.diff = X.diff[X.diff$t>=2005 & X.diff$t<=2040,]
  do.margin(X.diff,'part','abs') + facet_grid(cols=vars(case.lab)); fig.save(uid,'inf-diff-part-abs',w=10,h=3);
  do.margin(X.diff,'from','abs') + facet_grid(cols=vars(case.lab)); fig.save(uid,'inf-diff-from-abs',w=10,h=3);
  do.margin(X.diff,'to','abs')   + facet_grid(cols=vars(case.lab)); fig.save(uid,'inf-diff-to-abs',w=10,h=3);
}
# --------------------------------------------------------------------------------------------------
# main.ims.poster()
# main.foi()
# main.art()
