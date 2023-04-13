library('reshape2')
source('post/config.r')
source('post/wiw.r')

main.wiw = function(){
  X = clean.wiw.data(read.csvs('art-rf','wiw-diff','art',skip='base'))
  X = X[X$t > 2003 & X$t <= 2030,]
  for (margin in c('part','from','to')){
    g = do.margin(X,margin,type='abs',strat='case.lab') +
      facet_grid('~case.lab',scales='free_y') +
      labs(y='Additional Infections (\'000s)') +
      guides(colour=guide_legend(nrow=1),fill=guide_legend(nrow=1))
    fig.save(uid,N$sam,'art.wiw',margin, w=8,h=3)
  }
}

plot.obj.1 = function(X){
  out.labs = list(
    cuminfect = 'Cumulative Additional Infections (%)',
    incidence = 'Additional Incidence (%)')
  X$out = factor(X$out,names(out.labs),out.labs)
  g = plot.expo.box(expo.qs(X),unlist(out.labs),'case.lab',seq(2005,2030,5)) +
    facet_grid('out',scales='free_y') + 
    scale_color_manual(values=set.cols$art) +
    scale_fill_manual(values=set.cols$art) +
    labs(y='Relative Additional Infections (%)',x='Year',color='Left behind:',fill='Left behind:') +
    theme(legend.pos=c(.005,.99),legend.justification=c(0,1))
  fig.save(uid,N$sam,'art.1.rai',w=6,h=5)
}

num.obj.1 = function(X,t.hor=2030){
  X.hor = X[X$t==t.hor,]
  print(expo.qs(X.hor,q=c(.025,.5,.975)))
  x.case = function(case){ X.hor[X.hor$case==case,] }
  X.ref = x.case('fsw-cli-')
  X.ref$value = (X.ref$value - x.case('fsw+cli+')$value) / X.ref$value
  print(expo.qs(X.ref,q=c(.025,.5,.975)))
}

main.obj.1 = function(){
  X.wide = read.csvs('art-rf','expo','art',rdata='load')
  X = melt.expo.s(X.wide,pop='all',out=c('incidence','cuminfect'))
  base.value = X[X$case=='base',]$value
  X$value = 100 * (X$value - base.value) / base.value
  X = X[X$case!='base',]
  num.obj.1(X)
  plot.obj.1(X)
}

# main.wiw()
# main.obj.1()

