library('gifski')
source('utils/ops.r')
source('utils/plot.r')

N = 256
uid = '2022-01-24'
X = read.csv(root.path('data','mid',uid,paste0('sankey_base_N=0-',N-1,'.csv')))
X = X[X$infections>0,]

save_gif(lapply(unique(X$t),function(t){ # print(t)
  print(sankey(X[X$t==t,],norm=TRUE) + ggtitle(t))
}),'sankey.gif',width=600,height=800,delay=1/10,loop=TRUE)

# BUG: doesn't work as expected
# library('gganimate')
# g = sankey(X,norm=TRUE) + ggtitle('{frame_time}') + transition_time(t) + ease_aes('linear')
# animate(g,fps=10,w=600,h=800,renderer=gifski_renderer())
# anim_save('sankey.gif')
