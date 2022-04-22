source('post/config.r')

plot.posterior = function(X.long,vars=NULL){
  if (!is.null(vars)){ X.long = X.long[X.long$variable %in% vars,] }
  X.long$ll[!is.finite(X.long$ll)] = 2*min(X.long$ll,na.rm=TRUE)
  X.long$ll.cat = cut(X.long$ll,c(-Inf,quantile(X.long$ll,c(.99),na.rm=TRUE),Inf))
  # X.ll = X.long[X.long$variable=='t0_hiv',]; X.ll$variable = 'll'; X.ll$value = log(-X.ll$ll); # DEBUG
  # X.long = rbind(X.long,X.ll) # DEBUG
  g = ggplot(X.long,aes(x=value,fill=ll.cat,color=ll.cat)) +
    geom_density(alpha=.3) +
    # geom_histogram(bins=32) +
    scale_fill_viridis_d(option='inferno',begin=.2,end=.8) +
    scale_color_viridis_d(option='inferno',begin=.2,end=.8) +
    guides(fill='none',color='none') +
    labs(x='Parameter Value',y='Density') +
    theme_light() +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  if (length(vars) == 1){
    g = g + labs(x=vars) +
      theme(axis.title.x=element_text(family='mono'))
  } else {
    g = g + facet_wrap(vars(X.long$variable),scales='free') +
      theme(strip.text.x=element_text(family='mono'))
  }
}

# pre-compute once & save:
# X = read.csv(gen.name('cal','Ps','base',b='all'))
# X.long = melt(X,id=c('seed','ll'))
# save(X.long,file=gen.name('cal','Ps','case',b='all',ext='.rdata'))
load(file=gen.name('cal','Ps','case',b='all',ext='.rdata'))
n.var = length(unique(X.long$variable))
# plot
plot.posterior(X.long)
fig.save(uid,paste0('post-base-all'),w=2*sqrt(n.var),h=1.5*sqrt(n.var))
plot.posterior(X.long,vars='t0_hiv')
fig.save(uid,paste0('post-base-t0_hiv'),w=4,h=3)