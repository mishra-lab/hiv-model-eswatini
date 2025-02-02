library('reshape2')
library('geepack')
source('post/config.r')
source('post/wiw.r')

out.labs = list(
  cuminfect = 'Cumulative Additional Infections',
  incidence = 'Additional Incidence',
  diagnosed = 'Diagnosed\namong PLHIV',
  treated_c = 'Treated\namong Diagnosed',
  vls_c     = 'VLS\namong Treated',
  treated_u = 'Treated\namong PLHIV',
  vls_u     = 'VLS\namong PLHIV')

pop.labs = list(
  all = 'Overall',
  w   = 'Women Overall',
  m   = 'Men Overall',
  aq  = 'All Others',
  fsw = 'FSW',
  cli = 'Clients')

y.vars = c(
  'Cumulative Additional Infections' = 'cai',
  'Additional Incidence Rate' = 'air')
adj.vars = c(
  'DU Overall' = 'Du.all')
pred.vars = c(
  'd FSW'     = 'du.fsw',
  'd Clients' = 'du.cli')
mod.vars = c(
  'FSW % Population'     = 'px.fsw',
  'FSW / Women HIV IR'   = 'ir.fsw',
  'FSW Turnover'         = 'tur.fsw',
  'Clients % Population' = 'px.cli',
  'Clients / Men HIV IR' = 'ir.cli',
  'Clients Turnover'     = 'tur.cli')

var.lab.fun = function(vars){
  var.labs = c(y.vars,adj.vars,pred.vars,mod.vars)
  for (v in 1:length(var.labs)){
    vars = gsub(var.labs[v],names(var.labs[v]),vars)
  }
  vars = gsub('^.*:','',vars)
}

melt.expo.labs = function(X.wide,...,s=100){
  X = melt.expo.i(X.wide,...)
  X$value = s * X$value
  X$out = factor(X$out,names(out.labs),out.labs)
  X$pop = factor(X$pop,names(pop.labs),pop.labs)
  return(X)
}

melt.expo.cascade = function(X.wide,...){
  X = melt.expo.labs(X.wide,...,
    pop=c('all','aq','fsw','cli'),
    out=c('diagnosed','treated_c','vls_c','treated_u','vls_u'))
  return(X)
}

plot.clean.art = function(g){
  g = g + labs(x='Year',color='Left behind:',fill='Left behind:',linetype='Left behind:') +
    scale_color_manual(values=set.cols$art) +
    scale_fill_manual(values=set.cols$art) +
    scale_linetype_manual(values=set.lts$art)
  return(g)
}

plot.1.rai = function(X){
  X$out = factor(X$out,names(out.labs),out.labs)
  g = plot.expo.box(expo.qs(X),unlist(out.labs),'case.lab',seq(2005,2020,5)) +
    facet_grid('out',scales='free_y') +
    labs(y='Relative Additional Infections (%)')
  g = plot.clean.art(g)
  fig.save(uid,nid,'art.1.rai',w=6,h=5)
}

num.1.rai = function(X,t.hor=2020,s=1,rnd=1){
  X$value = s * X$value
  X.hor = X[X$t==t.hor,]
  print(expo.qs(X.hor,q=q3,trans=function(x){ round(x,rnd) }))
  x.case = function(case){ X.hor[X.hor$case==case,] }
  X.ref = x.case('fsw+cli+')
  X.ref$value = (x.case('fsw-cli-')$value - X.ref$value) / X.ref$value
  print(expo.qs(X.ref,q=q3,trans=function(x){ round(100*x,1) }))
}

main.1.rai = function(){
  X.wide = read.csvs('art-rf','expo','art',rdata='load')
  X = melt.expo.i(X.wide,pop='all',out=c('incidence','cuminfect'))
  base.value = X[X$case=='base',]$value
  num.1.rai(X[X$out=='cuminfect',],s=1,   rnd=0) # '000s
  num.1.rai(X[X$out=='incidence',],s=1000,rnd=0) # per 1000 PY
  X$value = (X$value - base.value) / base.value
  X = X[X$case!='base',]
  num.1.rai(X,s=100,rnd=1)
  plot.1.rai(X)
}

main.1.wiw = function(){
  X = clean.wiw.data(read.csvs('art-rf','wiw-diff','art',skip='base'))
  X = X[X$t > 2003 & X$t <= 2020,]
  for (margin in c('ptr','from','to')){
    g = wiw.margin(X,margin,type='abs',strat='case.lab') +
      facet_grid('~case.lab',scales='free_y') +
      labs(y='Additional Infections (\'000s)') +
      guides(colour=guide_legend(nrow=1),fill=guide_legend(nrow=1))
    fig.save(uid,nid,'art.1.wiw',margin,w=8,h=3)
  }
}

main.1.expo = function(){
  X.wide = read.csvs('art-rf','expo','art',rdata='load')
  X.wide = X.wide[X.wide$t >= 1995 & X.wide$t <= 2020,]
  X.cascade = melt.expo.cascade(X.wide)
  g = plot.expo.ribbon(expo.qs(X.cascade),unlist(out.labs),'case.lab') + facet_grid('out ~ pop') +
    lims(y=c(0,100)) + labs(y='Cascade Step (%)') + scale_x_continuous(breaks=c(2000,2010,2020))
  g = plot.clean.art(g); fig.save(uid,nid,'art.1.cascade',w=8,h=8)
  X.inc = melt.expo.labs(X.wide,out='incidence',pop='all',s=1000)
  g = plot.expo.ribbon(expo.qs(X.inc),unlist(out.labs),'case.lab',alpha=.1) +
    lims(y=c(0,NA)) + labs(y='HIV Incidence (per 1000 person-years)')
  g = plot.clean.art(g); fig.save(uid,nid,'art.1.inc',w=6,h=3)
}

main.2.cascade = function(){
  X.wide = read.csvs('art-ss','expo','sens',rdata='load')
  X = melt.expo.cascade(X.wide,t=2020)
  print(expo.qs(X[X$pop=='Overall',],q=q3,trans=function(x){ round(x,1) }))
  g = ggplot(expo.qs(X)) +
    geom_boxplot(q.aes('box','pop',x='pop'),stat='identity',alpha=.2,width=.8) +
    facet_grid('out') + coord_flip() +
    scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100)) +
    scale_color_manual(values=set.cols$pop.art) +
    scale_fill_manual(values=set.cols$pop.art) +
    labs(y='Cascade step (%)',x='Population',color='',fill='')
  g = plot.clean(g,legend.position='none')
  fig.save(uid,nid,'art.2.cascade',w=5,h=7)
}

load.2.data = function(rdata=''){
  X.list = list(
  'out' = do.call(rbind,lapply(c('base','sens'),function(case){
    X.case = filter.cols(read.csvs('art-ss','expo',case,rdata=rdata),
      out = c('prevalence','incidence','cuminfect','vls_u'))
    if (case=='base'){ X.case$ss = NA }
    return(X.case)
  })),
  'par' = read.csvs('art-ss','P0s','sens'))
}

clean.2.data = function(X.list,t.hor=2020,t.ref=2020,y.pop='all'){
  i.list = lapply(X.list,grep.i.col)
  x.out = function(...){ c(t(filter.cols(X.list$out,...)[,i.list$out])) }
  x.par = function(par){ c(t(filter.cols(X.list$par,par=par)[,i.list$par])) }
  X = data.frame(
    id      = colnames(X.list$par[i.list$par]),
    cai     = (x.out(case='sens',pop=y.pop,t=t.hor,out='cuminfect') - x.out(case='base',pop=y.pop,t=t.hor,out='cuminfect')) /
               x.out(case='base',pop=y.pop,t=t.hor,out='cuminfect'),
    air     = (x.out(case='sens',pop=y.pop,t=t.hor,out='incidence') - x.out(case='base',pop=y.pop,t=t.hor,out='incidence')) /
               x.out(case='base',pop=y.pop,t=t.hor,out='incidence'),
    Du.all  =  x.out(case='base',pop='all',t=t.ref,out='vls_u') - x.out(case='sens',pop='all',t=t.ref,out='vls_u'),
    du.fsw  =  x.out(case='sens',pop='all',t=t.ref,out='vls_u') - x.out(case='sens',pop='fsw',t=t.ref,out='vls_u'),
    du.cli  =  x.out(case='sens',pop='all',t=t.ref,out='vls_u') - x.out(case='sens',pop='cli',t=t.ref,out='vls_u'),
    ir.fsw  =  x.out(case='sens',pop='fsw',t=2010,out='incidence') / x.out(case='sens',pop='w',t=2010,out='incidence'),
    ir.cli  =  x.out(case='sens',pop='cli',t=2010,out='incidence') / x.out(case='sens',pop='m',t=2010,out='incidence'),
    tur.fsw = 1/x.par('dur_fsw'),
    tur.cli = 1/x.par('dur_cli'),
    px.fsw  = x.par('PX_fsw'),
    px.cli  = x.par('PX_cli'))
  X = X[order(X$id),]
}

fit.2.glm = function(X,y,std=TRUE){
  all.vars = c(adj.vars,pred.vars,mod.vars)
  if (std){ X[,all.vars] = apply(X[,all.vars],2,function(x){ (x-mean(x))/sd(x) }) }
  # if (std){ X[,all.vars] = apply(X[,all.vars],2,function(x){ (x-median(x))/iqr(x) }) } # DEBUG
  # if (std){ X[,all.vars] = apply(X[,all.vars],2,function(x){ rank(x)/length(x) }) } # DEBUG
  term.fun = function(...){ paste('(',paste(c(...),collapse=' + '),')') }
  f = paste(y,'* 100 ~ 1 +',
    term.fun(adj.vars),'+',
    term.fun(pred.vars),'+',
    term.fun(pred.vars),':',term.fun(mod.vars))
  # print(f) # DEBUG
  M = geeglm(formula(f),'gaussian',X,id=factor(X$id),corstr='e')
  # print(M); print(anova(M)) # DEBUG
  return(M)
}

plot.2.glm.effects = function(M.list,o.lab){
  f.lab = c('d Effects','d Interaction',paste(names(pred.vars),'Effect Modification'))
  E = do.call(rbind,lapply(names(M.list),function(name){
    e = summary(M.list[[name]])$coefficients
    e$var  = rownames(e)
    e$facet = factor(sapply(e$var,function(v){
      i = which(sapply(pred.vars,grepl,v))
      if      (length(i) == 0){ return(NA) }
      else if (length(i)  > 1){ return(f.lab[2]) }
      else if (!grepl(':',v)) { return(f.lab[1]) }
      else { return(f.lab[i+2]) }
    }),levels=f.lab)
    e$Est.low  = e$Estimate - e$Std.err * qnorm(.975)
    e$Est.high = e$Estimate + e$Std.err * qnorm(.975)
    e$var = var.lab.fun(e$var)
    e$name = name
    e = e[!is.na(e$facet),]
  }))
  dodge = position_dodge(width=.7)
  clrs = c(clr.map(6,'Blues')[c(5:3,6)],clr.map(6,'Reds')[c(6,5:3)]) # HACK
  g = ggplot(E,aes(y=var,x=100*Estimate,xmin=100*Est.low,xmax=100*Est.high,color=var)) +
    facet_grid('facet',scales='free',space='free') +
    geom_vline(xintercept=0,color=rgb(.8,.8,.8),lwd=1) +
    geom_point(position=dodge,size=1.5) +
    geom_errorbar(position=dodge,width=.25,lwd=.7) +
    scale_color_manual(values=clrs) +
    labs(x=paste('Effect on',o.lab,'(%) per SD change in variable     '),y='')
  g = plot.clean(g,legend.position='none')
}

plot.2.glm.resid = function(M,name){
  rfun = function(r){ r = abs(r); r = pmin(r,quantile(r,.99)) }
  X = cbind(M$data,data.frame(y=M$y,yp=predict(M),r=residuals(M)))
  g = ggplot(X,aes(x=y,y=yp,color=rfun(r))) +
    geom_abline(color='gray') +
    geom_point(shape=1,size=1) +
    scale_color_viridis(option='inferno',end=.9) +
    coord_fixed() +
    labs(x=paste('Simulation',name,'(%)'),
         y=paste('Regression',name,'(%)'),
         color='Residual')
  g = plot.clean(g,legend.position='top')
}

main.2.stats = function(){
  X.list = load.2.data('load')
  X = clean.2.data(X.list)
  # print(cor(as.matrix(X[,2:ncol(X)]))) # DEBUG
  for (name in names(y.vars)){
    y = y.vars[[name]]
    M = fit.2.glm(X,y)
    g = plot.2.glm.resid(M,name)
    fig.save(uid,nid,'art.2',y,'r',w=4,h=4)
    g = plot.2.glm.effects(list(x=M),name)
    fig.save(uid,nid,'art.2',y,w=6,h=5)
  }
}

# main.1.expo()
# main.1.rai()
# main.1.wiw()
# main.2.cascade()
# main.2.stats()
