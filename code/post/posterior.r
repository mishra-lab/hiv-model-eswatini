source('post/config.r')

q.cut = .01
id.vars = c('seed','ll','ll.cut')

pre.load = function(case){
  P = read.csv(gen.name('sam','Ps',case))
  P$ll[!is.finite(P$ll)] = -Inf
  P$foi_mode = NULL
  s.post = P$ll >= quantile(P$ll,1-q.cut,na.rm=TRUE)
  P.cor = cor(P[s.post,-which(names(P) %in% id.vars)],method='s')
  save(P,P.cor,file=gen.name('sam','Ps',case,b='all',ext='.rdata'))
}

cut.ll = function(P){
  labs = paste(c('Bottom','Top'),sprintf('%g%%',100*c(1-q.cut,q.cut)))
  ll.cut = cut.q(P$ll,1-q.cut,labels=labs)
}

plot.posterior = function(P,bins=32){
  vars = setdiff(names(P),id.vars)
  P$ll.cut = cut.ll(P)
  P.s = split(P,P$ll.cut)
  ksvp = unlist(lapply(vars,function(var){
    vp = paste(var,p.stars(ks.test(P.s[[1]][[var]],P.s[[2]][[var]])$p.value))
  }))
  P.long = melt(P,id.vars=id.vars)
  P.long$label = rep(ksvp,each=nrow(P))
  g = ggplot(P.long) +
    geom_freqpoly(aes(x=value,color=ll.cut,after_stat(density)),bins=bins) +
    facet_wrap('label',scales='free',ncol=7) +
    scale_color_viridis(option='inferno',discrete=TRUE,begin=0,end=.5) +
    scale_x_continuous(labels=function(x){ gsub('0\\.','.',x) }) +
    labs(x='Parameter Value',y='Density',color='Likelihood') +
    theme_light()
  g = plot.clean(g,legend.position='top',
    panel.grid.major = element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
}

plot.cor = function(X.cor,thr=.1,rect=0){
  i.thr = apply(abs(X.cor)-diag(nrow(X.cor)) > thr,1,any) # symmetric
  X.cor = X.cor[i.thr,i.thr]
  corrplot(X.cor,method='color',order='hclust',hclust.method='ward.D2',
    addrect=rect,rect.col='gray',
    col = scales::brewer_pal(palette='RdBu',direction=-1)(11),
    tl.cex=5/sqrt(nrow(X.cor)),tl.col='black')
}

# TODO: refactor with base.r etc.?
case = 'base'
# pre.load(case)
load(file=gen.name('sam','Ps',case,b='all',ext='.rdata')) # P, P.cor
g = plot.posterior(P); ggsave('post.distr.pdf',w=12,h=16)
plot.cor(P.cor,thr=.2,rect=7); file.rename('Rplots.pdf','post.cor.pdf')
