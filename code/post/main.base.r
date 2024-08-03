source('post/config.r')
source('post/post.r')
source('post/wiw.r')

main.post = function(){
  X = load.post.data(rdata='load')
  plot.post.cor(X$cor,thr=0)
  file.rename('Rplots.pdf','post.cor.pdf') # MAN
  g = plot.post.uni(X$pp,ncol=7) + scale_color_manual(values=c('#000000',clr))
  fig.save(uid,nid,'post.distr',w=12,h=16)
}

main.ll = function(){
  plot.ll.tform(); fig.save('','toy','ll.tf',w=4,h=6) # TOY
  X = read.csvs('imis','Ps','base',rdata='load')
  pid = read.csvs('fit','Ps','base',rdata='load')$id
  X$post = paste(X$batch,X$imis,X$id,sep='.') %in% pid
  X = rbind(
    cbind(X[X$imis==0,],g='Initial Samples'),
    cbind(X[X$imis >0,],g='IMIS Samples'),
    cbind(X[X$post   ,],g='Posterior Samples'))
  plot.ll.hist(X,y=after_stat(ndensity)) + facet_grid('g')
  fig.save(uid,nid,'ll.hist',w=5,h=4)
}

main.num = function(){
  X = melt.expo.i(read.csvs('art-ss','expo','base',rdata='load'))
  xo = filter.cols
  tf = function(f,d=1){ function(x){ round(x*f,d) } }
  xo.vs = function(X,pop1,pop2,...){
    X0 = xo(X,pop=pop1,...)
    X0$value = X0$value / xo(X,pop=pop2,...)$value
    X0$pop = paste0(pop1,':',pop2)
    return(X0) }
  print(expo.qs(xo(X,out='prevalence',pop='all',t=2020),q=q3,trans=tf(100)))
  print(expo.qs(xo(X,out='incidence', pop='all',t=2020),q=q3,trans=tf(1000)))
  print(expo.qs(xo.vs(X,'fsw','w',out='prevalence',t=2020),q=q3,trans=tf(1,2)))
  print(expo.qs(xo.vs(X,'cli','m',out='prevalence',t=2020),q=q3,trans=tf(1,2)))
}

main.wiw = function(){
  X = clean.wiw.data(read.csvs('fit','wiw','base',rdata='load'))
  X$infs.prop = X$infs / aggregate(infs~t,X,sum)$infs
  print(aggregate(infs.prop~ptr+t,X,sum)) # NUM
  wiw.ratio(X);         fig.save(uid,nid,'wiw.base.ratio',w=5,h=5)
  wiw.margin(X,'ptr');  fig.save(uid,nid,'wiw.base.ptr',  w=5,h=5.5)
  wiw.margin(X,'from'); fig.save(uid,nid,'wiw.base.from', w=5,h=8)
  wiw.margin(X,'to');   fig.save(uid,nid,'wiw.base.to',   w=5,h=8)
  wiw.alluvial.facet(X,seq(1990,2020,10),nrow=1); fig.save(uid,nid,'wiw.base.alluvial',w=10,h=5)
}

# main.num()
# main.post()
# main.ll()
# main.wiw()
