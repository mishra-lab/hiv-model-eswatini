suppressPackageStartupMessages({
  library('ggplot2')
  library('ggalluvial')
  library('viridis')
  library('reshape2')
})

fig.save = function(uid,...,w=7,h=7){
  fig.name = root.path('out','fig',uid,paste0(paste(...,sep='_'),'.pdf'))
  print(paste('saving:',fig.name))
  ggsave(fig.name,w=w,h=h)
}
