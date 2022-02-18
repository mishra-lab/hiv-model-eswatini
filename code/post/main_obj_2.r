library('geepack')
source('post/config.r')

clean.data = function(X.raw,t.hor=2040,t.cas=2020,t.prev=2000){
  X.base = X.raw[X.raw$case=='base',]
  X      = X.raw[X.raw$case!='base',]
  X$r.id = 1:N$rand
  X = X[order(X$r.id,X$seed),]
  o = function(Xi,output,t,pop){
    if (missing(t)){ t = t.hor }
    if (missing(pop)){ pop = 'all' }
    out = Xi[[paste0(output,'_',pop,'_',t)]]
  }
  X$PX_aq = 1 - X$PX_fsw - X$PX_cli
  X$dur_fsw = .8 * X$dur_fsw_l + .2 * X$dur_fsw_h
  X$prev_all  = o(X,'prev',t.prev,'all')
  X$prev_fsw  = o(X,'prev',t.prev,'fsw')
  X$prev_cli  = o(X,'prev',t.prev,'cli')
  X$prev_aq   = o(X,'prev',t.prev,'aq')
  X$pr_fsw.wq = o(X,'prev_ratio',t.prev,'fsw.wq')
  X$pr_cli.mq = o(X,'prev_ratio',t.prev,'cli.mq')
  X$ipr_all   = o(X,'inc',t.prev,'all') / X$prev_all
  X$inc.red = (o(X,'inc') - o(X.base,'inc')) / o(X,'inc')
  X$inc.add = (o(X,'inc') - o(X.base,'inc')) / o(X.base,'inc')
  X$inf.red = (o(X,'cuminf') - o(X.base,'cuminf')) / o(X,'cuminf')
  X$inf.add = (o(X,'cuminf') - o(X.base,'cuminf')) / o(X.base,'cuminf')
  X$d.vls.fsw = o(X.base,'vls_u',t.cas,'fsw') - o(X,'vls_u',t.cas,'fsw')
  X$d.vls.cli = o(X.base,'vls_u',t.cas,'cli') - o(X,'vls_u',t.cas,'cli')
  X$d.vls.aq  = o(X.base,'vls_u',t.cas,'aq')  - o(X,'vls_u',t.cas,'aq')
  # for (output in c('inc','cuminf')){ X[grepl(output,names(X))] = NULL } # remove raw outputs
  return(X)
}
plot.cascade = function(X){
  X.long = melt(X,m=grep('^diagnosed_|^treated_|^vls_',names(X)))
  X.long$variable = gsub('_c_','.c_',gsub('_u_','.u_',X.long$variable))
  X.long[,c('step','group')] = do.call(rbind,strsplit(as.character(X.long$variable),'_'))
  X.long$step = factor(X.long$step,levels=names(steps),labels=lab$steps)
  X.long$group = factor(X.long$group,levels=names(groups),labels=lab$group)
  # print(aggregate(value~group+step,X.long,quantile,c(.025,.5,.975)))
  g = ggplot(X.long,aes(x=100*value,fill=group,color=group)) + 
    geom_density(alpha=.2) +
    facet_grid(cols=vars(step)) +
    scale_color_manual(values=clr$groups) +
    scale_fill_manual(values=clr$groups) +
    labs(x='Cascade step value (%) in 2020',y='Density',color='Population',fill='Population') +
    theme_light() +
    theme(legend.position='top')
  return(g)
}
do.glm = function(X,y='inf.red',pred.vars,mod.vars,interact=TRUE,std=TRUE){
  pred.terms = iterms(pred.vars,ifelse(interact,2,1))
  mod.terms = paste('(',iterms(pred.vars,1),') : (',iterms(mod.vars,1),')')
  f = paste(y,'~ 0 +',pred.terms,'+',mod.terms)
  print(f)
  vars = c(pred.vars,mod.vars)
  X.m = X[,c(y,'seed',vars)]
  X.m$seed = factor(X.m$seed)
  if (std){ X.m[,vars] = apply(X.m[,vars],2,function(x){ (x-mean(x))/sd(x) }) }
  m = geeglm(formula(f),'gaussian',X.m,id=X.m$seed,corstr='i')
  # print(summary(m))
  return(m)
}
clean.var.labs = function(vars,v.labs){
  for (v in names(v.labs)){
    vars = gsub(v,paste0('(',v.labs[[v]],')'),vars)
  }
  vars = gsub(':',' x ',vars)
}
plot.effects = function(models,pred.vars,mod.vars,x.lab,y.lab='Variable',c.lab='Year'){
  w = qnorm(.975)
  m1 = length(models)==1
  f.lab = c('dVS Effects','dVS Interac.',paste(unlist(pred.vars),'Effect Mod'))
  E = do.call(rbind,lapply(names(models),function(name){
    e = summary(models[[name]])$coefficients
    e$vars = rownames(e)
    e$facet = sapply(e$vars,function(v){
      i = which(sapply(names(pred.vars),grepl,v))
      if (length(i) > 1){
        return(f.lab[2])
      } else if (!grepl(':',v)) {
        return(f.lab[1])
      } else {
        return(f.lab[i+2])
      }
    })
    e$facet = factor(e$facet,levels=f.lab)
    e$vars = clean.var.labs(e$vars,c(pred.vars,mod.vars))
    e$est.low  = e$Estimate - w*e$Std.err
    e$est.high = e$Estimate + w*e$Std.err
    e$name = name
    return(e)
  }))
  E$name = factor(E$name,levels=names(models),labels=names(models))
  g = ggplot(E,aes(y=vars,x=Estimate,xmin=est.low,xmax=est.high,color=name)) +
    geom_vline(xintercept=0,color=rgb(.8,.8,.8),lwd=1) +
    geom_point(size=1.4,position=position_dodge(width=.7)) +
    geom_linerange(lwd=.7,position=position_dodge(width=.7)) +
    scale_colour_viridis_d(option='inferno',begin=.2,end=.8,direction=ifelse(m1,1,-1)) +
    scale_y_discrete(limits=rev) +
    labs(x=x.lab,y=y.lab,color=c.lab) +
    facet_grid(rows=vars(facet),scales='free_y',space='free_y') +
    theme_light()
  if (m1){ g = g + theme(legend.position='none') }
  return(g)
}
effects.figsize = function(models){
  N.v = length(coef(models[[1]]))
  N.t = length(models)
  return(list(
    ts = N.t > 1,
    h = .5 + .04 * (6+N.t) * N.v,
    w = 7 + (N.t > 1)
  ))
  
}
plot.points = function(X,y,ylab,...){
  # TODO: facet by effects?
  X.long = melt(X,m=c('d.vls.fsw','d.vls.cli','d.vls.aq'))
  levels(X.long$variable) = c('dVS FSW','dVS Clients','dVS Lower Risk')
  g = ggplot(X.long,aes_string(x='100 * value',y=paste('100 *',y),...)) +
    geom_line(aes(group=seed),alpha=.3) +
    # geom_point(alpha=.5) +
    facet_grid(cols=vars(variable)) +
    scale_colour_viridis(option='inferno',begin=.1,end=.9) +
    labs(x='Reduction in Virally Suppressed (dVS)',y=ylab) +
    theme_light()
  return(g)
}
# config: this analysis
specs = list(
  'inc'   = list(y.def='inc.add',t=2040,y.lab='Additional Incidence'),
  'inf'   = list(y.def='inf.add',t=2040,y.lab='Cumulative Additional Infections'),
  'inf_t' = list(y.def='inf.add',t=c(2020,2030,2040),y.lab='Cumulative Additional Infections')
)
pred.vars = list(
  'd.vls.fsw' = 'dVS FSW',
  'd.vls.cli' = 'dVS Clients',
  'd.vls.aq'  = 'dVS Lower Risk')
mod.vars = list(
  'PX_fsw'      = 'FSW % Pop.',
  # 'prev_fsw'    = 'FSW HIV Prev.',
  # 'hiv_fsw'     = 'FSW % PLHIV',
  'dur_fsw_l'   = 'FSW Duration',
  'pr_fsw.wq'   = 'FSW / LR Women HIV PR',
  'PX_cli'      = 'Clients % Pop.',
  # 'prev_cli'    = 'Client HIV Prev.',
  # 'hiv_cli'     = 'Clients % PLHIV',
  'dur_cli'     = 'Client Duration',
  'pr_cli.mq'   = 'Clients / LR Men HIV PR'
  # 'A_swq_cli'   = 'Client Monthly FSW Visits',
  # 'prev_aq'     = 'Lower Risk HIV Prev.'
  # 'EHY_acute'   = 'HIV Acute EHY',
  # 'P_gud_fsw_l' = 'FSW LR GUD Prevlence',
  # 'A_mc'        = 'Main Yearly Sex Acts',
  # 'ipr_all'     = 'HIV IPR'
)

X.raw = load.csvs('sens',case.rand(),batches=1:10)
X = clean.data(X.raw)
g = plot.cascade(X); fig.save(uid,'obj_2_cascade',w=10,h=3)

for (s in names(specs)){
  spec = specs[[s]]
  models = lapply(setNames(spec$t,spec$t),function(t){
    do.glm(clean.data(X.raw,t),spec$y.def,names(pred.vars),names(mod.vars))
  })
  g = plot.effects(models,pred.vars=pred.vars,mod.vars=mod.vars,
    x.lab=paste('Effect on',spec$y.lab),y.lab='Standardized dVS Variable')
  fs = effects.figsize(models)
  if (!fs$ts){ print(summary(models[[1]])) }
  fig.save(uid,'obj_2',s,w=fs$w,h=fs$h)
}
# g = plot.points(X,y.def,y.lab,color='ipr_all'); ggsave('Rplots.pdf',w=12,h=4); q() # DEBUG
# q()