library('gifski')
source('utils/ops.r')
source('utils/plot.r')

# TODO: rename waifw / alluvial & update dates + data

N = 1024
uid = '2022-02-02'
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
# --------------------------------------------------------------------------------------------------
load.data = function(){
  X = read.csv(root.path('data','mid',uid,paste0('sankey_base_N=0-',N-1,'.csv')))
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
  ggsave('alluvial.pdf',w=12,h=16)
}
do.alluvial.gif = function(X){
  save_gif(lapply(unique(X$t),function(t){ # print(t)
    print(alluvial(X[X$t==t,],norm=TRUE) + ggtitle(t))
  }),'alluvial.gif',width=600,height=800,delay=1/10,loop=TRUE)
}
do.margin = function(X,margin){
  # partnership types
  X.abs = aggregate(formula(paste('infections ~ t +',margin)),X,sum)
  X.rel = do.norm(X.abs,sum(X.abs$t==X.abs$t[1]))
  X.abs$scale = 'Absolute (\'000s)'
  X.rel$scale = 'Proportion (%)'
    if (margin == 'part'){ print(X.rel[X.rel$t %in% seq(1990,2040,5),]) }
  if (margin %in% c('from','to')){
    i.men = grepl('Men\\s|\\sClients',X.abs[[margin]])
    X.abs$infections[i.men] = -X.abs$infections[i.men]
  }
  X.long = rbind(X.abs,X.rel)
  g = ggplot(X.long,aes_string(x='t',y='infections',color=margin,fill=margin)) +
    geom_area(alpha=.7,size=.1) +
    facet_wrap(vars(scale),scales='free_y') +
    scale_color_manual(values=clrs[[margin]]) +
    scale_fill_manual(values=clrs[[margin]]) +
    labs(x='Year',y='Infections',color='',fill='') +
    theme_light()
  if (margin %in% c('from','to')){ g = g + scale_y_continuous(labels=abs) }
  return(g)
}

# --------------------------------------------------------------------------------------------------
X = load.data()
do.margin(X,'part'); ggsave('inf-part.pdf',w=8,h=4)
do.margin(X,'to');   ggsave('inf-to.pdf'  ,w=8,h=4)
do.margin(X,'from'); ggsave('inf-from.pdf',w=8,h=4)
do.alluvial.facet(X)
do.alluvial.gif(X)

# BUG: doesn't work as expected
# library('gganimate')
# g = alluvial(X,norm=TRUE) + ggtitle('{frame_time}') + transition_time(t) + ease_aes('linear')
# animate(g,fps=10,w=600,h=800,renderer=gifski_renderer())
# anim_save('alluvial.gif')
