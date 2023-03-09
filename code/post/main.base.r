source('post/config.r')
# N$sam = 1000 # DEBUG

base.wiw = function(){
  source('post/wiw.r')
  X = clean.wiw.data(read.csvs('fit','wiw','base','all'))
  do.ratio(X);         fig.save(uid,N$sam,'wiw.base.ratio',w=5,h=5);
  do.margin(X,'part'); fig.save(uid,N$sam,'wiw.base.part',w=5,h=5.5);
  do.margin(X,'from'); fig.save(uid,N$sam,'wiw.base.from',w=5,h=8);
  do.margin(X,'to');   fig.save(uid,N$sam,'wiw.base.to',  w=5,h=8);
  do.alluvial.facet(X,seq(1990,2040,10),nrow=2); fig.save(uid,N$sam,'wiw.base.alluvial',w=8,h=10)
}

# base.wiw()
