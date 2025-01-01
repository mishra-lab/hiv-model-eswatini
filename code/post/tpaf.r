plot.tpaf = function(X,t.hor,t='t.hor',grp='case.lab',pop='all',mode='box',...){
  Xe = filter.cols(X,t.hor=t.hor,pop=pop,...)
  Xe$pop = factor(Xe$pop,levels=names(strat.labs),labels=strat.labs)
  Xe$tpaf.path = factor(Xe$tpaf.path,levels=names(strat.labs),labels=strat.labs)
  args = list(X=Xe,out='cuminfect',grp=grp,scale=100)
  if (mode=='box'){    plot.fun = plot.expo.box;    args$X$t = Xe[[t]]; args$t = unique(Xe[[t]]) }
  if (mode=='ribbon'){ plot.fun = plot.expo.ribbon; args$X$t = Xe[[t]] }
  g = do.call(plot.fun,args) + lims(y=c(0,NA)) + ylab('TPAF (%)')
}
