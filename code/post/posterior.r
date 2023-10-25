source('post/config.r')

id.vars = c('id','batch','imis','ll','lp','wt','case','case.lab','distr')
i.vars = function(X){ intersect(colnames(X),id.vars) }
p.vars = function(X){ setdiff(colnames(X),id.vars) }

load.post.data = function(case='base',rdata='load'){
  X.imis = read.csvs('imis','Ps',case,b=seq(N$batch)-1,rdata=rdata)
  X.post = read.csvs('fit', 'Ps',case,rdata=rdata)
  X = list(
    pp = rbind(
      cbind(X.imis[X.imis$imis==0,],distr='Prior'),
      cbind(X.post,distr='Posterior')),
    cor = cor(X.post[,p.vars(X.post)],method='p'))
}

plot.ll.bi = function(X,ll.min=-1e6){
  Xll = aggregate(ll~batch+imis,X,mean)
  Xll = rbind(
    cbind(Xll,g='Log Likelihood', value=rescale(pmax(Xll$ll,ll.min))),
    cbind(Xll,g='Rank Likelihood',value=rescale(rank(Xll$ll))))
  g = ggplot(Xll,aes(x=factor(imis),y=factor(batch),fill=value,color=value)) +
    facet_grid('g') +
    geom_tile() +
    scale_x_discrete(labels=NULL,breaks=NULL) +
    scale_y_discrete(labels=NULL,breaks=NULL) +
    scale_fill_viridis(option='inferno') +
    scale_color_viridis(option='inferno') +
    labs(x='IMIS Iteration',y='IMIS Batch')
  g = plot.clean(g,legend.position='none')
}

plot.ll.hist = function(X,ll.min=-1e3){
  X$ll[X$ll < ll.min] = NA
  hargs = list(binwidth=10,color='white',lwd=.1,alpha=.5)
  g = ggplot(map=aes(x=ll,y=after_stat(ndensity))) +
    do.call(geom_histogram,c(hargs,list(data=X[X$post==FALSE,]))) +
    do.call(geom_histogram,c(hargs,list(data=X[X$post==TRUE,],fill=clr))) +
    scale_y_continuous(labels=NULL) +
    labs(x='Log Likelihood',y='Normalized Density')
  g = plot.clean(g)
}

plot.post.uni = function(X,color='distr',ncol=NULL,bins=32,stats=FALSE,p.thr=1){
  X.long = melt(X,id.vars=i.vars(X)) 
  X.long$variable = factor(X.long$variable,levels=sort(levels(X.long$variable)))
  if (stats){
    p = ad.tests(X,color)
    X.long = merge(X.long,cbind(p=p,p.star=sapply(p,p.stars),variable=p.vars(X)))
    X.long$variable = paste(X.long$variable,X.long$p.star)
    X.long = X.long[X.long$p <= p.thr,]
  }
  g = ggplot(X.long,aes_string(color=color)) +
    geom_freqpoly(aes(x=value,after_stat(density)),bins=bins) +
    stat_summary(aes(x=value,y=0),fun=mean,geom='point',shape=1,orientation='y') +
    facet_wrap('variable',scales='free',ncol=ncol) +
    scale_x_continuous(labels=function(x){ gsub('0\\.','.',x) }) +
    labs(x='Parameter Value',y='Density',color='')
  g = plot.clean(g,legend.position='top',
    panel.grid.major=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
}

plot.post.cor = function(X.cor,thr=.1){
  i.thr = apply(abs(X.cor)-diag(nrow(X.cor)) > thr,1,any) # symmetric
  X.cor = X.cor[i.thr,i.thr]
  corrplot(X.cor,method='color',order='hclust',hclust.method='ward.D2',
    col = scales::brewer_pal(palette='RdBu',direction=-1)(11),
    tl.cex=30/nrow(X.cor),tl.col='black')
}

ad.tests = function(X,g='distr'){
  p.vals = par.lapply(p.vars(X),function(var){ print(var)
    p.val = kSamples::ad.test(split(X[[var]],X[[g]]))$ad[1,3]
  })
}
