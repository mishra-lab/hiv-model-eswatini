source('post/config.r')
source('post/wiw.r')
# N$sam = 1000 # DEBUG

base.ll = function(){
  X = read.csvs('fit','ll','base')
  g = plot.hist(X,x='ll',bw=10,color='case',fill='case') +
    scale_color_manual(values='white') +
    scale_fill_manual(values='#CC0033') +
    labs(x='Log Likelihood',y='Count') +
    theme(legend.position='none')
  fig.save(uid,N$sam,'fit.ll.base',w=5,h=3)
}

base.wiw = function(){
  X = clean.wiw.data(read.csvs('fit','wiw','base'))
  do.ratio(X);         fig.save(uid,N$sam,'wiw.base.ratio',w=5,h=5);
  do.margin(X,'part'); fig.save(uid,N$sam,'wiw.base.part', w=5,h=5.5);
  do.margin(X,'from'); fig.save(uid,N$sam,'wiw.base.from', w=5,h=8);
  do.margin(X,'to');   fig.save(uid,N$sam,'wiw.base.to',   w=5,h=8);
  do.alluvial.facet(X,seq(1990,2040,10),nrow=2); fig.save(uid,N$sam,'wiw.base.alluvial',w=8,h=10)
}

# base.wiw()
# base.ll()
