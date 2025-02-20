source('post/config.r')

serv.labs  = list(base='None',px='Prevention',tx='Treatment',ax='Both')
tline.labs = list(base='Base',opt='Optimistic',pes='Pessimistic',non='No Recovery')
serv.clrs = c(base='#000000',px='#CC0066',tx='#00CCCC',ax='#FFCC00')

out.scales = list(
  incidence  = 1000,
  prevalence = 100,
  cuminfect  = 1000,
  cumdeath   = 1000,
  diagnosed  = 100,
  treated_c  = 100,
  vls_c      = 100,
  treated_u  = 100,
  vls_u      = 100)

zero.cum = function(X.wide,out,t0=2025){
  i = grep.i.col(X.wide)
  bo = X.wide$out == out
  bt = X.wide$t >= t0
  X.wide[bo & !bt,i] = 0
  X.ref = X.wide[bo & X.wide$t==t0,i]
  for (t in unique(X.wide$t[bt])){
    j = which(bo & X.wide$t==t)
    X.wide[j,i] = X.wide[j,i] - X.ref
  }
  return(X.wide)
}

add.vs = function(X.all,sub,fun,op='vs'){
  X.ref = subset(X.all,eval(sub))
  X.ref = rbind.lapply(1:(nrow(X.all)/nrow(X.ref)),function(i){ X.ref })
  i = grep.i.col(X.all)
  X.all.vs = X.all
  X.all.vs[,i] = fun(X.all[,i],X.ref[,i])
  X.all.vs$op = op
  X.all = rbind(X.all,X.all.vs)
}

clean.future = function(X.wide,t0=2020,ti=2025,tf=2035){
  X.wide = subset(X.wide, t >= t0 & t <= tf)
  X.wide = zero.cum(X.wide,'cuminfect',t0=ti)
  X.wide = zero.cum(X.wide,'cumdeath', t0=ti)
  X.wide$serv  = factor(X.wide$serv, names(serv.labs), serv.labs)
  X.wide$tline = factor(X.wide$tline,names(tline.labs),tline.labs)
  i = grep.i.col(X.wide)
  for (out in names(out.scales)){
    j = which(X.wide$out == out)
    X.wide[j,i] = X.wide[j,i] * out.scales[[out]]
  }
  X.wide = add.vs(X.wide,
    quote(serv==serv.labs$base),
    function(all,ref){ all-ref })
}

plot.future.box = function(X,t.ref=2035){
  out.labs = list(cuminfect='Infections',cumdeath='Deaths')
  X = melt.expo.i(X.wide,pop='all',op='vs',t=t.ref,
    out=names(out.labs),serv=serv.labs[2:4])
  X$out = factor(X$out,names(out.labs),out.labs)
  g = ggplot(X,aes(x='',y=value,color=serv)) +
    facet_grid('out~tline',scales='free_y') +
    scale_y_continuous(trans='identity') +
    scale_color_manual(values=unname(serv.clrs),drop=FALSE) +
    labs(x='',y='Additional HIV Outcomes 2025 - 2035',color='Disruption') +
    geom_boxplot(fill=NA,outlier.shape=1,width=1)
  g = plot.clean(g)
}

plot.future.inc = function(X.wide,pop='all'){
  X = melt.expo.i(X.wide,pop=pop,op='raw',out='incidence')
  i = which(X$serv==serv.labs$base)
  X.base = rbind(X[i,],X[i,],X[i,])
  X.base$tline = rep(unlist(tline.labs[2:4]),each=nrow(X.base)/3)
  g = ggplot(X[-i,],aes(x=t,y=value,color=serv)) +
    facet_grid('~tline') + ylim(c(0,NA)) +
    scale_color_manual(values=unname(serv.clrs)) +
    labs(x='Year',y='HIV Incidence (per 1000 person-years)',color='Disruption') +
    stat_summary(data=X.base,fun='median',geom='line') +
    stat_summary(fun='median',geom='line')
  g = plot.clean(g)
}

plot.future.vls = function(X.wide,t.ref=2035){
  pops = c('Overall'='all','Female Sex Workers'='fsw')
  X = melt.expo.i(X.wide,out='vls_u',op='raw',pop=pops,t=t.ref)
  X$pop = factor(X$pop,pops,names(pops))
  i = lapply(names(pops),function(pop){ which(X$pop==pop) })
  X$value = 100 - X$value
  X$value[i[[1]]] = -X$value[i[[1]]]
  g = ggplot(X,aes(x=value,y=serv,fill=serv)) +
    facet_grid('tline~pop',scales='free',space='free') +
    scale_fill_manual(values=unname(serv.clrs)) +
    scale_y_discrete(limits=rev,drop=TRUE) +
    scale_x_continuous(labels=abs) +
    stat_summary(fun='median',geom='bar',width=.8) +
    stat_summary(fun.min=min,fun.max=max,geom='errorbar',color='#999999',width=.4) +
    labs(x='% PLHIV Not Virally Suppressed in 2035',y='',fill='Disruption')
  g = plot.clean(g)
}


sets$future = set.labs$future = c('base','txo','pxo','axo','txp','pxp','axp','txn','pxn','axn')
X.wide = clean.future(read.csvs('future','expo','future'))
g = plot.future.inc(X.wide,pop='all'); fig.save(uid,nid,'future.inc.all',w=8,h=3)
g = plot.future.inc(X.wide,pop='fsw'); fig.save(uid,nid,'future.inc.fsw',w=8,h=3)
g = plot.future.vls(X.wide);           fig.save(uid,nid,'future.vls.bar',w=6,h=4)
g = plot.future.box(X.wide);           fig.save(uid,nid,'future.box',w=6,h=4)
