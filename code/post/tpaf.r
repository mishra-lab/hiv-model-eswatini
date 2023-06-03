
plot.tpaf = function(X,t.hor,grp='case.lab',pop='all',mode='box',...){
  Xe = filter.cols(X,t.hor=t.hor,pop=pop,...)
  Xe$pop = factor(Xe$pop,levels=names(slice.labs),labels=slice.labs)
  Xe$tpaf.pop = factor(Xe$tpaf.pop,levels=names(slice.labs),labels=slice.labs)
  args = list(X=Xe,out='cuminfect',grp=grp,scale=100)
  if (mode=='box'){ plot.fun = plot.expo.box; args$X$t = args$X$t.hor; args$t = t.hor }
  if (mode=='ribbon'){ plot.fun = plot.expo.box }
  g = do.call(plot.fun,args) + lims(y=c(0,NA))
}

plot.tpaf.box = function(X,t.hor=c(1,3,10),...){
  g = plot.tpaf(X,t.hor=t.hor,...) +
    labs(x='Time Horizon (years)',y='TPAF (%)')
  g = plot.clean(g,legend.position='top')
}
