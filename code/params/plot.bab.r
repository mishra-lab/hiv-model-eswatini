library('ggplot2')
library('viridis')
eps = .001
x = seq(eps,1-eps,eps)
N = c(10,20,40)
p = c(.01,.1,.50)
data = do.call(rbind,lapply(N,function(Ni){
  do.call(rbind,lapply(p,function(pi){
    cbind(N=Ni,p=pi,rbind(
      data.frame(x=(0:Ni)/Ni,d=dbinom(0:Ni,Ni,pi),distr='Binom'),
      data.frame(x=x,d=dbeta(x,Ni*pi,Ni*(1-pi))/Ni,distr='Beta')
    ))
}))}))
g = ggplot(mapping=aes(x=x,y=d,color=distr)) +
  geom_line(data=data[data$distr=='Beta',]) +
  geom_step(data=data[data$distr=='Binom',],direction='mid') +
  facet_grid(p~N,scales='free_y',labeller=label_both) +
  scale_color_viridis(discrete=TRUE,option=3,begin=0,end=.5) +
  scale_x_continuous(breaks=c(.25,.5,.75)) +
  labs(x='p',y='Density',color='') +
  theme_light() +
  theme(legend.position='top',axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave('Rplots.pdf',w=6,h=6)
