source('post/config.r')

id.vars = c('id','batch','imis','foi_mode','ll','lp','wt','case','case.lab','distr')
i.vars = function(X){ intersect(colnames(X),id.vars) }
p.vars = function(X){ setdiff(colnames(X),id.vars) }

load.post.data = function(case='base',rdata='load'){
  X.imis = read.csvs('imis','Ps',case,b=seq(N$batch)-1,rdata=rdata)
  X.post = read.csvs('fit', 'Ps',case,rdata=rdata)
  X = list(
    pp = rbind(
      cbind(X.imis[X.imis$imis==0,],distr='Prior'),
      cbind(X.post,distr='Posterior')),
    cor = cor(X.post[,p.vars(X.post)],method='s'))
}

plot.ll.hist = function(X,bw=5,ll.min=-1e3,...){
  X$ll[X$ll < ll.min] = NA
  g = ggplot(X,aes(x=ll,...)) +
    geom_histogram(binwidth=bw,lwd=.1,color='white') +
    lims(x=c(ll.min,0)) +
    labs(x='Log Likelihood',y='Density')
  g = plot.clean(g)
}

plot.post.uni = function(X,color='distr',ncol=NULL,bins=32,stats=FALSE,p.thr=1){
  X.long = melt(X,id.vars=i.vars(X)) 
  X.long$variable = factor(X.long$variable,levels=sort(levels(X.long$variable)))
  if (stats){
    p = distr.tests(X,color)
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

plot.post.cor = function(X.cor,thr=0){
  i.thr = apply((abs(X.cor)-diag(nrow(X.cor))) > thr,1,any) # symmetric
  corrplot(X.cor[i.thr,i.thr],
    type='upper',method='color',diag=FALSE,
    col = scales::brewer_pal(palette='RdBu',direction=-1)(11),
    tl.cex=3.33/sqrt(nrow(X.cor)),tl.col='black')
}

distr.tests = function(X,group='distr'){
  p.vals = par.lapply(p.vars(X),function(var){ pout(var)
    p.val = kSamples::qn.test(split(X[[var]],X[[group]]))$qn[2]
    # p.val = kSamples::ad.test(split(X[[var]],X[[group]]))$ad[1,3] # DEBUG
  })
}

plot.ll.tform = function(){
  lls = -1*10^seq(1,6,.01)
  tform = function(x,q=.1){ quantile(x,1-q)/x }
  X = rbind(
    data.frame(ll=lls,y=-log10(-lls),g='Log Log Likelihood'),
    data.frame(ll=lls,y=lls/1000,g='Log Likelihood / 1000'),
    data.frame(ll=lls,y=1e3*exp(lls),g='Likelihood x 1000'),
    data.frame(ll=lls,y=tform(lls),g='Transform'))
  g = ggplot(X,aes(x=-log10(-ll),y=y)) +
    geom_line(color=clr) +
    geom_line(data=data.frame(ll=quantile(lls,c(0,.9,.9)),y=c(1,1,0),g='Transform'),lty=2) +
    facet_grid('g',scales='free') +
    scale_x_continuous(breaks=seq(-6,-1),labels=paste0('-10^',seq(6,1))) +
    labs(x='Log Likelihood',y='Value')
  plot.clean(g,axis.text.x=ggtext::element_markdown())
}
