source('post/config.r')
source('post/post.r')
source('post/wiw.r')

main.post = function(){
  X = load.post.data()
  plot.post.cor(X$cor,thr=0)
  file.rename('Rplots.pdf','post.cor.pdf')
  g = plot.post.uni(X$pp,ncol=7) + scale_color_manual(values=c('#000000',clr))
  ggsave('post.distr.pdf',w=12,h=16) # MAN
}

main.ll = function(){
  plot.ll.tform(); fig.save('','toy','ll.tf',w=4,h=6) # DEBUG
  X = read.csvs('imis','Ps','base',rdata='load')
  pid = read.csvs('fit','Ps','base',rdata='load')$id
  X$post = paste(X$batch,X$imis,X$id,sep='.') %in% pid
  plot.ll.hist(X); fig.save(uid,nid,'ll.hist',w=5,h=4)
}

main.wiw = function(){
  X = clean.wiw.data(read.csvs('fit','wiw','base'))
  do.ratio(X);         fig.save(uid,nid,'wiw.base.ratio',w=5,h=5)
  do.margin(X,'part'); fig.save(uid,nid,'wiw.base.part', w=5,h=5.5)
  do.margin(X,'from'); fig.save(uid,nid,'wiw.base.from', w=5,h=8)
  do.margin(X,'to');   fig.save(uid,nid,'wiw.base.to',   w=5,h=8)
  do.alluvial.facet(X,seq(1990,2040,10),nrow=2); fig.save(uid,nid,'wiw.base.alluvial',w=8,h=10)
}

# main.post()
# main.ll()
# main.wiw()
