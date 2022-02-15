library('gifski')
source('utils/ops.r')
source('utils/plot.r')

# TODO: rename waifw / alluvial & update dates + data

N = 1024
uid = '2022-02-13'
pops = list(
  'Women Low'  = list(clr=rgb(1.0,.60,.60)),
  'Women Med'  = list(clr=rgb(.80,.40,.40)),
  'LR FSW'     = list(clr=rgb(.60,.20,.20)),
  'HR FSW'     = list(clr=rgb(.40,.00,.00)),
  'Men Low'    = list(clr=rgb(.60,.60,1.0)),
  'Men Med'    = list(clr=rgb(.40,.40,.80)),
  'LR Clients' = list(clr=rgb(.20,.20,.60)),
  'HR Clients' = list(clr=rgb(.00,.00,.40)))
pop.labs = names(pops)
parts = list(
  'Main / Spousal' = list(clr=rgb(.26,.04,.41)),
  'Casual'         = list(clr=rgb(.58,.15,.40)),
  'New Sex Work'   = list(clr=rgb(.87,.32,.23)),
  'Reg Sex Work'   = list(clr=rgb(.99,.65,.04)))
part.labs = names(parts)
clrs = list(
  from = unname(sapply(pops, function(pop){pop$clr})),
  to   = unname(sapply(pops, function(pop){pop$clr})),
  part = unname(sapply(parts,function(part){part$clr}))
)
case.defs = list(
  'FSW-Cli+'=list(clr=rgb(.9,.3,.3),title='Left Behind: FSW'),
  'FSW+Cli-'=list(clr=rgb(.3,.3,.9),title='Left Behind: Clients'),
  'FSW-Cli-'=list(clr=rgb(.8,.3,.8),title='Left Behind: FSW & Clients'),
  'FSW+Cli+'=list(clr=rgb(.9,.7,.0),title='Left Behind: Neither'))
cases = names(case.defs)
case.titles = unname(sapply(case.defs,function(case){case$title}))
# --------------------------------------------------------------------------------------------------
load.data = function(diff=FALSE){
  if (diff){
    which = 'infs-diff'
  } else {
    which = 'infs'
    cases = c(cases,'base')
    case.titles = c(case.titles,'Base Case')
  }
  X = do.call(rbind,lapply(cases,function(case){
    X.i = read.csv(root.path('data','mid',uid,paste0(which,'_',case,'_N=0-',N-1,'.csv')))
    X.i$case = case
    return(X.i)
  }))
  X$title = factor(X$case,levels=cases,labels=case.titles)
  X = X[X$infections>0,] # for speed
  # paste hack to maintain order in alluvial
  X$from = factor(1+X$fs*4+X$fi,levels=1:8,labels=pop.labs)
  X$to   = factor(1+X$ts*4+X$ti,levels=1:8,labels=pop.labs)
  X$part = factor(1+X$p,levels=1:4,labels=part.labs)
  return(X)
}
do.norm = function(X,n,scale=100){
  total.inf = aggregate(infections~t,X,sum)$infections
  X$infections = scale * X$infections / rep(total.inf,n)
  return(X)
}
alluvial = function(X,fill='part',norm=FALSE){
  X$from = factor(1+X$fs*4+X$fi,levels=1:8,labels=paste(pop.labs,''))
  X$to   = factor(1+(1-X$ts)*4+X$ti,levels=1:8,labels=paste0('',c(pop.labs[5:8],pop.labs[1:4])))
  X.long = to_lodes_form(X,key='pop',axis=c('from','to'))
  if (norm){ X.long = do.norm(X.long,sum(X$t==X$t[1]),200) }
  ylab = paste('Yearly Infections',ifelse(norm,'(%)','(\'000s)'))
  g = ggplot(X.long,aes(y=infections,x=pop,stratum=stratum,label=stratum,alluvium=alluvium)) +
    geom_alluvium(aes_string(fill=fill),width=.4) +
    geom_stratum(color=rgb(.5,.5,.5),alpha=.5,width=.4) +
    geom_text(stat='stratum',size=3) +
    labs(x='Population',y=ylab,fill='') +
    scale_x_discrete(expand=c(.15,.15)) +
    scale_fill_manual(values=clrs[[fill]]) +
    theme_light()
  return(g)
}
do.alluvial.facet = function(X){
  tvec = seq(1990,2040,10)
  g = alluvial(X[X$t %in% tvec,],norm=TRUE) +
    facet_wrap(~factor(t)) +
    theme(legend.position='top')
}
do.alluvial.gif = function(X){
  save_gif(lapply(unique(X$t),function(t){ # print(t)
    print(alluvial(X[X$t==t,],norm=TRUE) + ggtitle(t))
  }),'alluvial.gif',width=600,height=800,delay=1/10,loop=TRUE)
}
do.margin = function(X,margin,rel=TRUE){
  # partnership types
  X.abs = aggregate(formula(paste('infections ~ t + case + title +',margin)),X,sum)
  X.rel = do.norm(X.abs,sum(X.abs$t==X.abs$t[1]))
  X.abs$scale = 'Absolute (\'000s)'
  X.rel$scale = 'Proportion (%)'
  if (margin == 'part'){ print(X.rel[X.rel$t %in% seq(1990,2040,5),]) }
  if (margin %in% c('from','to')){
    i.men = grepl('Men\\s|\\sClients',X.abs[[margin]])
    X.abs$infections[i.men] = -X.abs$infections[i.men]
  }
  if (rel){ X.long = rbind(X.abs,X.rel) } else { X.long = X.abs }
  g = ggplot(X.long,aes_string(x='t',y='infections',color=margin,fill=margin)) +
    geom_area(alpha=.7,size=.1) +
    scale_color_manual(values=clrs[[margin]]) +
    scale_fill_manual(values=clrs[[margin]]) +
    labs(x='Year',y='Infections',color='',fill='') +
    theme_light()
  if (margin %in% c('from','to')){ g = g + scale_y_continuous(labels=abs) }
  if (rel){ g = g + facet_wrap(vars(scale),scales='free_y') }
  return(g)
}

# --------------------------------------------------------------------------------------------------
X = load.data()
X.base = X[X$case=='base',]
do.margin(X.base,'part');  fig.save(uid,'inf-part',w=8,h=4)
do.margin(X.base,'to');    fig.save(uid,'inf-to'  ,w=8,h=4)
do.margin(X.base,'from');  fig.save(uid,'inf-from',w=8,h=4)
do.alluvial.facet(X.base); fig.save(uid,'inf-alluvial',w=12,h=16)
do.alluvial.gif(X.base)
X.diff = load.data(diff=TRUE)
do.margin(X.diff,'part',rel=FALSE) + facet_grid(cols=vars(title)); fig.save(uid,'inf-diff-part',w=10,h=3)
do.margin(X.diff,'to',  rel=FALSE) + facet_grid(cols=vars(title)); fig.save(uid,'inf-diff-to',  w=10,h=3)
do.margin(X.diff,'from',rel=FALSE) + facet_grid(cols=vars(title)); fig.save(uid,'inf-diff-from',w=10,h=3)

# BUG: doesn't work as expected
# library('gganimate')
# g = alluvial(X,norm=TRUE) + ggtitle('{frame_time}') + transition_time(t) + ease_aes('linear')
# animate(g,fps=10,w=600,h=800,renderer=gifski_renderer())
# anim_save('alluvial.gif')
