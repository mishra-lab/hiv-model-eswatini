source('post/config.r')

plot.posterior = function(X.long){
  X.long$ll[!is.finite(X.long$ll)] = 2*min(X.long$ll,na.rm=TRUE)
  X.long$ll.cat = cut(X.long$ll,c(-Inf,quantile(X.long$ll,c(.99),na.rm=TRUE),Inf))
  # X.ll = X.long[X.long$variable=='t0_hiv',]; X.ll$variable = 'll'; X.ll$value = log(-X.ll$ll); # DEBUG
  # X.long = rbind(X.long,X.ll) # DEBUG
  g = ggplot(X.long,aes(x=value,fill=ll.cat,color=ll.cat)) +
    facet_wrap(vars(X.long$variable),scales='free') +
    geom_density(alpha=.3) +
    # geom_histogram(bins=32) +
    scale_fill_viridis_d(option='inferno',begin=.2,end=.8) +
    scale_color_viridis_d(option='inferno',begin=.2,end=.8) +
    guides(fill='none',color='none') +
    labs(x='Parameter Value',y='Density') +
    theme_light() +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
      strip.text.x=element_text(family='mono'))
}

# pre-compute once & save:
# X = read.csv(gen.name('cal','Ps','base',b='all'))
# X.long = melt(X,id=c('seed','ll'))
# save(X.long,file=gen.name('cal','Ps','case',b='all',ext='.rdata'))
load(file=gen.name('cal','Ps','case',b='all',ext='.rdata'))
# plot
plot.posterior(X.long[seq(1,nrow(X.long),1),])
ggsave('Rplots.pdf',w=18,h=12)