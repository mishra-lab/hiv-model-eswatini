source('post/config.r')

q.cut = .01
id.vars = c('seed','ll','distr','case','case.lab')
i.vars = function(X){ intersect(colnames(X),id.vars) }
p.vars = function(X){ setdiff(colnames(X),id.vars) }
s.post = function(X){ X$ll > quantile(X$ll,1-q.cut,na.rm=TRUE) }

pre.load = function(case){
  X = read.csv(gen.name('sam','Ps',case))
  X$ll[!is.finite(X$ll)] = -Inf
  X$foi_mode = NULL
  X.cor = cor(X[s.post(X),-which(names(X) %in% id.vars)],method='s')
  save(X,X.cor,file=gen.name('sam','Ps',case,ext='.rdata'))
}

def.pp = function(X){
  X = rbind(cbind(X,distr='Prior'),cbind(X[s.post(X),],distr='Posterior'))
}

plot.post.uni = function(X,color='distr',ncol=NULL,bins=32,stats=TRUE,p.thr=1){
  X.long = melt(X,id.vars=i.vars(X))
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

plot.post.cor = function(X.cor,thr=.1,rect=1){
  i.thr = apply(abs(X.cor)-diag(nrow(X.cor)) > thr,1,any) # symmetric
  X.cor = X.cor[i.thr,i.thr]
  corrplot(X.cor,method='color',order='hclust',hclust.method='ward.D2',
    addrect=rect,rect.col='gray',
    col = scales::brewer_pal(palette='RdBu',direction=-1)(11),
    tl.cex=5/sqrt(nrow(X.cor)),tl.col='black')
}

ad.tests = function(X,g='distr'){
  p.vals = par.lapply(p.vars(X),function(var){
    p.val = kSamples::ad.test(split(X[[var]],X[[g]]))$ad[1,3]
  })
}
