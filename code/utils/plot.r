suppressPackageStartupMessages({
  library('ggplot2')
  library('ggalluvial')
  library('corrplot')
  library('viridis')
  library('reshape2')
  library('scales')
})
clr = '#CC0033'
clr.map = RColorBrewer::brewer.pal
q.aes = function(q,grp,x='t',scale=1,...){
  qq = function(qi,op){ paste0(scale,'*',op,'(q',1-(1-qi)/2,',q',(1-qi)/2,')') }
  if (q==0){ return(aes_string(x=x,y=qq(0,'pmin'),color=grp,linetype=grp,...)) }
  if (is.numeric(q)){ return(aes_string(x=x,ymin=qq(q,'pmin'),ymax=qq(q,'pmax'),fill=grp,...)) }
  if (q=='box'){ return(aes_string(x=paste0('factor(',x,')'),...,
    ymin=qq(.95,'pmin'),ymax=qq(.95,'pmax'),lower=qq(.5,'pmin'),upper=qq(.5,'pmax'),
    middle=qq(0,'pmin'),color=grp,fill=grp,group=paste0('interaction(',x,',',grp,')'))) }
}
plot.expo.ribbon = function(X,out,grp,q=.9,alpha=.25,...){
  Xe = filter.cols(X,out=out)
  g = ggplot(Xe,aes(x=t)) +
    geom_ribbon(q.aes(q,grp=grp,...),alpha=alpha) +
    geom_line(q.aes(0,grp=grp,...))
  g = plot.clean(g)
}
plot.expo.box = function(X,out,grp,t,w=.75,...){
  Xe = filter.cols(X,out=out,t=t)
  g = ggplot(Xe,aes(x=factor(t))) +
    geom_boxplot(q.aes('box',grp=grp,...),stat='identity',
      alpha=.25,width=.8*w,position=position_dodge(width=w))
  g = plot.clean(g)
}
plot.hist = function(X,x,bw=1,...){
  g = ggplot(X,aes_string(x=x,...)) +
    geom_histogram(binwidth=bw,alpha=.6,lwd=.1)
  g = plot.clean(g)
}
fig.save = function(uid,nid,...,w=7,h=7){
  fig.name = root.path('out','fig',uid,nid,paste0(paste(...,sep='.'),'.pdf'))
  print(paste('saving:',fig.name))
  ggsave(fig.name,w=w,h=h)
}
plot.clean = function(g,...){
  g = g + theme_light() + theme(...,
    strip.background=element_rect(fill='gray85'),
    strip.text.x=element_text(color='black'),
    strip.text.y=element_text(color='black'))
}
