suppressPackageStartupMessages({
  library('ggalluvial')
  library('viridis')
})

sankey = function(X,fill='part',norm=FALSE){
  pop.labs  = c('W Low','W Med','LR FSW','HR FSW','M Low','M Med','LR Cli','HR Cli')
  part.labs = c('Main / Spousal','Casual','New Sex Work','Reg Sex Work')
  # paste hack to maintain order
  X$from = factor(1+X$fs*4+X$fi,levels=1:8,labels=paste(pop.labs,''))
  X$to   = factor(1+(1-X$ts)*4+X$ti,levels=1:8,labels=paste0('',c(pop.labs[5:8],pop.labs[1:4])))
  X$part = factor(1+X$p,levels=1:4,labels=part.labs)
  X.long = to_lodes_form(X,key='pop',axis=c('from','to'))
  if (norm){ X.long$infections = 200 * X.long$infections / sum(X.long$infections) }
  ylab = paste('Yearly Infections',ifelse(norm,'(%)','(\'000s)'))
  g = ggplot(X.long,aes(y=infections,x=pop,stratum=stratum,label=stratum,alluvium=alluvium)) +
    geom_alluvium(aes_string(fill=fill)) +
    geom_stratum(color='black',alpha=.5) +
    geom_text(stat='stratum') +
    labs(x='Population',y=ylab,fill='') +
    scale_x_discrete(expand=c(.15,.15)) +
    scale_fill_viridis_d(option='inferno',begin=.2,end=.8) +
    theme_light()
  return(g)
}