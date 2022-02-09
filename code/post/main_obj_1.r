library('ggplot2')
source('utils/ops.r')
source('utils/plot.r')
# --------------------------------------------------------------------------------------------------
# config: model fits
N    = 1024
uid  = '2022-02-02'
# config: cases
cases = list(# NOTE: order affects plot
  'FSW-Cli-AQ+'=list(clr=rgb(.7,.3,.7),short='--+',title='Leave Behind: FSW & Clients'),
  'FSW-Cli+AQ+'=list(clr=rgb(.9,.3,.3),short='-++',title='Leave Behind: FSW'),
  'FSW+Cli-AQ+'=list(clr=rgb(.3,.3,.9),short='+-+',title='Leave Behind: Clients'),
  'FSW-Cli-AQ-'=list(clr=rgb(.9,.6,.9),short='---',title='Leave Behind: Everybody'),
  'FSW-Cli+AQ-'=list(clr=rgb(.9,.6,.6),short='-+-',title='Leave Behind: FSW & Lower Risk'),
  'FSW+Cli-AQ-'=list(clr=rgb(.6,.6,.9),short='+--',title='Leave Behind: Clients & Lower Risk'),
  'FSW+Cli+AQ-'=list(clr=rgb(.8,.8,.8),short='++-',title='Leave Behind: Lower Risk'),
  'base'       =list(clr=rgb(.0,.0,.0),short='+++',title='base'))
colors = unname(sapply(cases,function(case){case$clr}))
titles = unname(sapply(cases,function(case){case$title}))
shorts = unname(sapply(cases,function(case){case$short}))
t.hors = seq(2000,2050,5)
t.labs = t.hors; t.labs[seq(2,length(t.labs),2)] = '' # HACK
# --------------------------------------------------------------------------------------------------
# load & process data
load.clean.data = function(){
  X.raw = do.call(rbind,lapply(names(cases),function(case){
    X.raw.i = read.csv(root.path('data','mid',uid,paste0('sens_',case,'_N=0-',N-1,'.csv')))
    X.raw.i$label = cases[[case]]$label
    return(X.raw.i)
  }))
  vars = c('cuminf','incid')
  X = reshape(X.raw,idvar=c('seed','case'),direction='long',times=t.hors,timevar='year',
    v.names=vars,varying=lapply(vars,function(i){ grepl(i,names(X.raw)) }))
  X = X[order(X$case,X$seed),]
  X.ref = X[X$case=='base',]
  X$inf.red = (X$cuminf - X.ref$cuminf) / X$cuminf
  X$inf.ext = (X$cuminf - X.ref$cuminf) / X.ref$cuminf
  X$inc.red = (X$incid  - X.ref$incid)  / X$incid
  X$inc.ext = (X$incid  - X.ref$incid)  / X.ref$incid
  X$title = X$case; levels(X$title) = titles
  X$short = X$case; levels(X$short) = shorts
  return(X)
}
plot.obj.1 = function(X,y='100*inf.red',ylab='Infections averted (%)',fill='title',...){
  g = ggplot(X[X$case!='base',],aes_string(x='factor(year)',y=y,fill=fill,...)) +
    geom_boxplot(show.legend=FALSE,outlier.size=.5,lwd=.3,width=.6) +
    facet_wrap(vars(title),ncol=3,dir='v') +
    labs(x='Year',y=ylab) +
    scale_x_discrete(labels=t.labs) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    theme_light()
  return(g)
}
save.tex = function(value,...,dec=3,per=1){
  value = value * per
  dec = dec - log10(per)
  str = sprintf(paste0('%.',dec,'f'),value)
  writeLines(str,root.path('out','tex',uid,...,create=TRUE))
}
save.med.ci = function(values,fname,dec=3,per=1){
  v = quantile(values,c(.025, .5, .975))
  save.tex(v[2],paste0(fname,''),   dec=dec,per=per)
  save.tex(v[1],paste0(fname,'.lo'),dec=dec,per=per)
  save.tex(v[3],paste0(fname,'.hi'),dec=dec,per=per)
}
numeric.obj.1 = function(X){
  X. = function(short,col){ X[X$short==short,][[col]] }
  Nf = length(unique(X.('+++','seed')))
  save.tex(N,'n.sample',dec=0)
  save.tex(10,'top.pct.fit',dec=0)
  save.tex(Nf,'n.fit',dec=0)
  save.tex(10,'n.rand',dec=0) # obj.2
  save.med.ci(X.('+++','prev_all'),'2020/prev.all',per=100)
  save.med.ci(X.('+++','inc_all'),'2020/inc.all',per=1000)
  save.med.ci(X.('+++','prev_ratio_fsw.wq'),'2020/pr.fsw')
  save.med.ci(X.('+++','prev_ratio_cli.mq'),'2020/pr.cli')
}
# --------------------------------------------------------------------------------------------------
X = load.clean.data()
# print(names(X))
# numeric.obj.1(X)
# X.t = X[X$year==2040,]
# f = function(i){X.t[X.t$title==titles[i],]$inc.red}
# hist((f(2) + f(3) + 0*f(7))/f(1))
# q()
plot.obj.1(X,y='100*inf.red',ylab='Infections averted (%) in base case vs ...');  fig.save(uid,'obj_1_inf_red')
plot.obj.1(X,y='100*inf.ext',ylab='Additional infections (%) vs base case');      fig.save(uid,'obj_1_inf_ext')
plot.obj.1(X,y='100*inc.red',ylab='Incidence reduction (%) in base case vs ...'); fig.save(uid,'obj_1_inc_red')
plot.obj.1(X,y='100*inc.ext',ylab='Additional incidence (%) vs base case');       fig.save(uid,'obj_1_inc_ext')