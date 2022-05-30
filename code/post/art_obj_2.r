library('geepack')
source('post/config.r')

clean.data = function(X.base,X,t.hor=2040,t.cas=2020,t.prev=2000){
  X.sens$r.id = 1:N$sens
  X.sens = X.sens[order(X.sens$r.id,X.sens$seed),]
  # print(all(X.sens$seed==X.base$seed)) # DEBUG
  o = function(Xi,output,pop,t=t.hor){
    return(Xi[[paste0(output,'_',pop,'_',t)]])
  }
  X$PX_aq = 1 - X$PX_fsw - X$PX_cli
  X$dur_fsw = .8 * X$dur_fsw_l + .2 * X$dur_fsw_h
  X$prev_all  = o(X,'prev','all',t.prev)
  X$prev_fsw  = o(X,'prev','fsw',t.prev)
  X$prev_cli  = o(X,'prev','cli',t.prev)
  X$prev_aq   = o(X,'prev','aq',t.prev)
  X$prev_wq   = o(X,'prev','wq',t.prev)
  X$prev_mq   = o(X,'prev','mq',t.prev)
  X$hiv_fsw   = o(X,'prev','fsw',t.prev) * X$PX_fsw
  X$hiv_cli   = o(X,'prev','cli',t.prev) * X$PX_cli
  X$hiv_aq    = o(X,'prev','aq',t.prev)  * X$PX_aq
  X$pr_fsw.wq = o(X,'prev','fsw',t.prev) / o(X,'prev','wq',t.prev)
  X$ir_fsw.wq = o(X,'inc','fsw',t.prev)  / o(X,'inc','wq',t.prev)
  X$pr_cli.mq = o(X,'prev','cli',t.prev) / o(X,'prev','mq',t.prev)
  X$ir_cli.mq = o(X,'inc','cli',t.prev)  / o(X,'inc','mq',t.prev)
  X$ipr_all   = o(X,'inc','all',t.prev)  / o(X,'prev','all',t.prev)
  for (pop in c('all','fsw','cli','aq','wq','mq')){
    X[[paste0('inc.red.',pop)]] = (o(X,'inc',pop) - o(X.base,'inc',pop)) / o(X,'inc',pop)
    X[[paste0('inc.add.',pop)]] = (o(X,'inc',pop) - o(X.base,'inc',pop)) / o(X.base,'inc',pop)
    X[[paste0('inf.red.',pop)]] = (o(X,'cuminf',pop) - o(X.base,'cuminf',pop)) / o(X,'cuminf',pop)
    X[[paste0('inf.add.',pop)]] = (o(X,'cuminf',pop) - o(X.base,'cuminf',pop)) / o(X.base,'cuminf',pop)
  }
  X$d.vls.all = o(X.base,'vls_u','all',t.cas) - o(X,'vls_u','all',t.cas)
  X$d.vls.fsw = o(X.base,'vls_u','fsw',t.cas) - o(X,'vls_u','fsw',t.cas)
  X$d.vls.cli = o(X.base,'vls_u','cli',t.cas) - o(X,'vls_u','cli',t.cas)
  X$d.vls.aq  = o(X.base,'vls_u','aq',t.cas)  - o(X,'vls_u','aq',t.cas)
  X$d.vls.wq  = o(X.base,'vls_u','wq',t.cas)  - o(X,'vls_u','wq',t.cas)
  X$d.vls.mq  = o(X.base,'vls_u','mq',t.cas)  - o(X,'vls_u','mq',t.cas)
  X$ru.vls.fsw = (1-o(X,'vls_u','fsw',t.cas)) / (1-o(X,'vls_u','all',t.cas))
  X$ru.vls.cli = (1-o(X,'vls_u','cli',t.cas)) / (1-o(X,'vls_u','all',t.cas))
  X$ru.vls.aq  = (1-o(X,'vls_u','aq',t.cas))  / (1-o(X,'vls_u','all',t.cas))
  X$ru.vls.wq  = (1-o(X,'vls_u','wq',t.cas))  / (1-o(X,'vls_u','all',t.cas))
  X$ru.vls.mq  = (1-o(X,'vls_u','mq',t.cas))  / (1-o(X,'vls_u','all',t.cas))
  return(X)
}
plot.cascade = function(X){
  X.long = melt(X,m=grep('^diagnosed_|^treated_|^vls_',names(X)))
  X.long$variable = gsub('_c_','.c_',gsub('_u_','.u_',X.long$variable))
  X.long[,c('step','group')] = do.call(rbind,strsplit(as.character(X.long$variable),'_'))
  X.long$step = factor(X.long$step,levels=names(spec$cascade),labels=sget('cascade','lab'))
  X.long$group = factor(X.long$group,levels=names(spec$groups),labels=sget('groups','lab'))
  print(aggregate(value~group+step,X.long,quantile,c(.025,.5,.975)))
  g = ggplot(X.long,aes(x=100*value,fill=group,color=group)) + 
    geom_density(alpha=.2) +
    facet_wrap(vars(step),ncol=5,scales='free') +
    scale_color_manual(values=sget('groups','clr'),na.translate=FALSE) +
    scale_fill_manual(values=sget('groups','clr'),na.translate=FALSE) +
    labs(x='Cascade step value (%) in 2020',y='Density',color='Population',fill='Population') +
    theme_light() +
    theme(legend.position='top',axis.ticks.y=element_blank(),axis.text.y=element_blank())
  return(g)
}
do.glm = function(X,y,adj.vars,pred.vars,mod.vars,std=TRUE){
  mod.terms = paste('(',iterms(mod.vars,1),') : (',iterms(pred.vars,1),')')
  f = paste(y,'~',iterms(adj.vars,1),'+',mod.terms,'+',iterms(pred.vars,1))
  vars = c(adj.vars,pred.vars,mod.vars)
  X.m = X[,c(y,'seed',vars)]
  X.m$seed = factor(X.m$seed)
  if (std){ X.m[,vars] = apply(X.m[,vars],2,function(x){ (x-mean(x))/sd(x) }) }
  m = geeglm(formula(f),'gaussian',X.m,id=X.m$seed,corstr='e')
  return(m)
}
clean.var.labs = function(vars,v.labs){
  for (v in names(v.labs)){
    vars = gsub(v,paste0('(',v.labs[[v]],')'),vars)
  }
  vars = gsub(':',' x ',vars)
}
plot.effects = function(models,pred.vars,mod.vars,x.lab,y.lab='Variable',vary=' '){
  c.lab = list(' '='.',t.hors='Year',groups='Among')[[vary]]
  w = qnorm(.975)
  f.lab = c('RU Effects','RU Interac.',paste(unlist(pred.vars),'Effect Mod'))
  E = do.call(rbind,lapply(names(models),function(name){
    e = summary(models[[name]])$coefficients
    e$vars = rownames(e)
    e$facet = sapply(e$vars,function(v){
      i = which(sapply(names(pred.vars),grepl,v))
      if (length(i) == 0){
        return(NA)
      } else if (length(i) > 1){
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
    e = e[!is.na(e$facet),]
    return(e)
  }))
  E$name = factor(E$name,levels=names(models),labels=names(models))
  g = ggplot(E,aes(y=vars,x=Estimate,xmin=est.low,xmax=est.high,color=name)) +
    geom_vline(xintercept=0,color=rgb(.8,.8,.8),lwd=1) +
    geom_point(size=1.4,position=position_dodge(width=.7)) +
    geom_linerange(lwd=.7,position=position_dodge(width=.7)) +
    scale_y_discrete(limits=rev) +
    labs(x=x.lab,y=y.lab,color=c.lab) +
    facet_grid(rows=vars(facet),scales='free_y',space='free_y') +
    theme_light()
  if (length(models)==1){
    g = g + theme(legend.position='none') +
      scale_color_manual(values=rgb(.69,.20,.35))
  } else {
    g = g + scale_color_manual(values=sget(vary,'clr'))
  }
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
# --------------------------------------------------------------------------------------------------
ipops = paste0('inf.add.',names(spec$groups)); names(ipops) = sget('groups','lab')
glm.specs = list(
  'inc'     = list(y.def='inc.add.all',t.hor=2030,y.lab='Additional Incidence'),
  'inf'     = list(y.def='inf.add.all',t.hor=2040,y.lab='Cumulative Additional Infections'),
  'inf_t'   = list(y.def='inf.add.all',t.hor=c(2020,2030,2040),y.lab='Cumulative Additional Infections',c.lab='Year'),
  'inf_pop' = list(y.def=ipops,t.hor=2040,y.lab='Cumulative Additional Infections',c.lab='Population')
)
adj.vars = list(
  'd.vls.all' = 'dU Overall'
)
pred.vars = list(
  'ru.vls.fsw' = 'RU FSW',
  'ru.vls.cli' = 'RU Clients'
)
mod.vars = list(
  'PX_fsw'      = 'FSW % Pop.',
  'pr_fsw.wq'   = 'FSW / LR Women HIV PR',
  'dur_fsw_l'   = 'FSW Duration',
  'PX_cli'      = 'Clients % Pop.',
  'pr_cli.mq'   = 'Clients / LR Men HIV PR',
  'dur_cli'     = 'Clients Duration'
)
# --------------------------------------------------------------------------------------------------
X.sens = read.big.csv(gen.name('art','keyout','sens','all'));
X.base = read.big.csv(gen.name('art','keyout','base','all'));
X = clean.data(X.base,X.sens)
g = plot.cascade(X.sens); fig.save(uid,'art_2_cascade',w=10,h=3)
for (s in names(glm.specs)){
  gs = glm.specs[[s]]
  if (length(gs$t.hor)>1){
    vary = 't.hors'
    models = lapply(setNames(gs$t.hor,gs$t.hor),function(t.hor){
      do.glm(clean.data(X.base,X.sens,t.hor),gs$y.def,
        names(adj.vars),names(pred.vars),names(mod.vars))
    })
  } else if (length(gs$y.def)>1){
    vary = 'groups'
    models = lapply(gs$y.def,function(y.def){
      do.glm(clean.data(X.base,X.sens,gs$t.hor),y.def,
        names(adj.vars),names(pred.vars),names(mod.vars))
    })
  } else {
    vary = ' '
    models = list('.'=do.glm(clean.data(X.base,X.sens,gs$t.hor),gs$y.def,
      names(adj.vars),names(pred.vars),names(mod.vars)))
  }
  g = plot.effects(models,pred.vars=pred.vars,mod.vars=mod.vars,vary=vary,
    x.lab=paste('Effect on',gs$y.lab),y.lab='Standardized Variable')
  fs = effects.figsize(models)
  if (!fs$ts){ print(summary(models[[1]])) }
  fig.save(uid,'art_2',s,w=fs$w,h=fs$h)
}
