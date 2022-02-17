library('geepack')
source('post/config.r')

clean.data = function(X.raw,t.hor=2040){
  X.base = X.raw[X.raw$case=='base',]
  X      = X.raw[X.raw$case!='base',]
  X$r.id = 1:N$rand
  X = X[order(X$r.id,X$seed),]
  o = function(Xi,output){
    out = Xi[[paste0(output,'_',t.hor)]]
  }
  X$inc.red = (o(X,'incid')  - o(X.base,'incid'))  / o(X,'incid')
  X$inc.add = (o(X,'incid')  - o(X.base,'incid'))  / o(X.base,'incid')
  X$inf.red = (o(X,'cuminf') - o(X.base,'cuminf')) / o(X,'cuminf')
  X$inf.add = (o(X,'cuminf') - o(X.base,'cuminf')) / o(X.base,'cuminf')
  X$d.vls.fsw = X.base$vls_u_FSW - X$vls_u_FSW
  X$d.vls.cli = X.base$vls_u_Cli - X$vls_u_Cli
  X$d.vls.aq  = X.base$vls_u_AQ  - X$vls_u_AQ
  X$PX_aq = 1 - X$PX_fsw - X$PX_cli
  X$hiv_fsw = X$PX_fsw * X$prev_fsw
  X$hiv_cli = X$PX_cli * X$prev_cli
  X$hiv_aq  = X$PX_aq * X$prev_aq
  X$dur_fsw = .8 * X$dur_fsw_l + .2 * X$dur_fsw_h
  X$d.vls.fsw.pp = X$d.vls.fsw / (X$prev_fsw * X$PX_fsw)
  X$d.vls.cli.pp = X$d.vls.cli / (X$prev_cli * X$PX_cli)
  X$d.vls.aq.pp  = X$d.vls.aq  / (X$prev_aq  * X$PX_aq)
  X$ipr_all   = X$inc_all / X$prev_all
  for (output in c('incid','cuminf')){ X[grepl(output,names(X))] = NULL } # remove raw outputs
  return(X)
}
plot.cascade = function(X){
  X.long = melt(X,m=grep('^diagnosed_|^treated_|^vls_',names(X)))
  X.long$variable = gsub('_c','.c',gsub('_u','.u',X.long$variable))
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
  if (std){
    vars = c(pred.vars,mod.vars)
    X[,vars] = apply(X[,vars],2,function(x){ (x-mean(x))/sd(x) })
  }
  print(f)
  m = geeglm(formula(f),'gaussian',X,id=factor(X$seed),corstr='i')
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
  f.lab = c('dVS Effects','dVS Interac.',paste(unlist(pred.vars),'Effect Mod.'))
  E = do.call(rbind,lapply(names(models),function(name){
    e = summary(models[[name]])$coefficients
    e$vars = rownames(e)
    e$facet = sapply(e$vars,function(v){
      # TODO: pop out independent factors from modifications
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
  g = ggplot(E,aes(y=vars,x=Estimate,xmin=est.low,xmax=est.high,color=name)) +
    geom_vline(xintercept=0,color=rgb(.8,.8,.8),lwd=1) +
    geom_point(size=1.4,position=position_dodge(width=.7)) +
    geom_linerange(lwd=.7,position=position_dodge(width=.7)) +
    scale_colour_viridis_d(option='inferno',begin=.2,end=.8) +
    scale_x_continuous(trans='pseudo_log') +
    scale_y_discrete(limits=rev) +
    labs(x=x.lab,y=y.lab,color=c.lab) +
    facet_grid(rows=vars(facet),scales='free_y',space='free_y') +
    theme_light()
  if (length(M) == 1){ g = g + theme(legend.position='none') }
  return(g)
}
effects.figsize = function(M){
  N.v = length(coef(M[[1]]))
  N.t = length(M)
  return(list(
    ts = N.t > 1,
    h = .5 + .04 * (6+N.t) * N.v,
    w = 7 + (N.t > 1)
  ))
  
}
plot.points = function(X,y,ylab,...){
  # TODO: facet by effects?
  X.long = melt(X,m=c('d.vls.fsw','d.vls.cli','d.vls.aq'))
  levels(X.long$variable) = c('NVS FSW','NVS Clients','NVS Lower Risk')
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
y.def = 'inf.add'
y.lab = 'Cumulative Additional Infections (%)'
t.def = 2040
pred.vars = list(
  'd.vls.fsw' = 'dVS FSW',
  'd.vls.cli' = 'dVS Clients',
  'd.vls.aq'  = 'dVS Lower Risk')
mod.vars = list(
  'PX_fsw'            = 'FSW % Pop.',
  # 'hiv_fsw'           = 'FSW % PLHIV',
  'dur_fsw_l'         = 'FSW Duration',
  'prev_ratio_fsw.wq' = 'FSW / Women HIV PR',
  'PX_cli'            = 'Clients % Pop.',
  # 'hiv_cli'           = 'Clients % PLHIV',
  'dur_cli'           = 'Client Duration',
  'prev_ratio_cli.mq' = 'Clients / Men HIV PR',
  # 'A_swq_cli'         = 'Client Monthly FSW Visits',
  # 'EHY_acute'         = 'HIV Acute EHY',
  # 'P_gud_fsw_l'       = 'FSW LR GUD Prevlence',
  # 'A_mc'              = 'Main Yearly Sex Acts',
  'ipr_all'           = 'HIV IPR'
)
X.raw = load.csvs('sens',case.rand())
X = clean.data(X.raw,t.def)
g = plot.cascade(X); fig.save(uid,'obj_2_cascade',w=10,h=3)
# g = plot.points(X,y.def,y.lab,color='ipr_all'); ggsave('Rplots.pdf',w=12,h=4); q() # DEBUG

for (t.vec in list(c(2020,2030,2040),2040)){
  M = lapply(setNames(t.vec,t.vec),function(t){
    do.glm(clean.data(X.raw,t),y.def,names(pred.vars),names(mod.vars))
  })
  g = plot.effects(M,pred.vars=pred.vars,mod.vars=mod.vars,
    x.lab=paste('Effect on',y.lab),y.lab='Standardized dVS Variable')
  fs = effects.figsize(M)
  if (!fs$ts){ print(summary(M[[1]])) }
  fig.save(uid,'obj_2',paste0('eff',ifelse(fs$ts,'_t','')),w=fs$w,h=fs$h)
}
