source('post/config.r')
source('post/posterior.r')
source('post/wiw.r')
# N$sam = 1000 # DEBUG

main.ll = function(){
  X = read.csvs('fit','ll','base')
  g = plot.hist(X,x='ll',bw=10,color='case',fill='case') +
    scale_color_manual(values='white') +
    scale_fill_manual(values='#CC0033') +
    labs(x='Log Likelihood',y='Count') +
    theme(legend.position='none')
  fig.save(uid,N$sam,'fit.ll.base',w=5,h=3)
}

main.post = function(){
  load(gen.name('sam','Ps','base',ext='.rdata'))
  plot.post.cor(X.cor,thr=.2,rect=7) # MAN
  file.rename('Rplots.pdf','post.cor.pdf')
  g = plot.post.uni(def.pp(X),ncol=7) + scale_color_manual(values=c('#000000','#CC0033'))
  ggsave('post.distr.pdf',w=12,h=16) # MAN
}

main.wiw = function(){
  X = clean.wiw.data(read.csvs('fit','wiw','base'))
  do.ratio(X);         fig.save(uid,N$sam,'wiw.base.ratio',w=5,h=5);
  do.margin(X,'part'); fig.save(uid,N$sam,'wiw.base.part', w=5,h=5.5);
  do.margin(X,'from'); fig.save(uid,N$sam,'wiw.base.from', w=5,h=8);
  do.margin(X,'to');   fig.save(uid,N$sam,'wiw.base.to',   w=5,h=8);
  do.alluvial.facet(X,seq(1990,2040,10),nrow=2); fig.save(uid,N$sam,'wiw.base.alluvial',w=8,h=10)
}

# main.ll()
# main.post()
# main.wiw()
