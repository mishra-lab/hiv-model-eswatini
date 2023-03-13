source('post/config.r')
source('post/wiw.r')
source('post/tpaf.r')
# N$sam = 1000 # DEBUG

foi.lab = 'FOI Approach'

plot.clean.foi = function(g,leg='top'){
  g = g + labs(color=foi.lab,fill=foi.lab,linetype=foi.lab) +
    scale_fill_manual(values=set.cols$foi) +
    scale_color_manual(values=set.cols$foi) +
    scale_linetype_manual(values=set.lts$foi) +
    theme(legend.position=leg)
}

plot.ep = function(Xepp,op,ylab,leg='top',oname='incidence'){
    g = plot.expo.ribbon(filter.cols(Xepp,op=op),oname,grp='case.lab') +
      facet_grid('pop ~ case.lab.x',scales='free_y') +
      labs(x='Year',y=ylab)
    g = plot.clean.foi(g,leg=leg)
}

main.ep = function(){
  # load & clean data
  X = read.csvs('foi-ep','expo','foi')
  Xe = filter.cols(X,pop=c('aq','fsw','cli'))
  Xe$pop = factor(Xe$pop,levels=names(slice.labs),labels=slice.labs)
  Xe$case.lab.x = Xe$case.lab
  # duplicate base data for each other case
  b.cols = grep('^q.*|case.lab$|^op$',names(Xe))
  b.rows = Xe$case=='base'
  Xep = Xe[!b.rows,]
  Xep[b.cols] = Xe[b.rows,b.cols]
  Xepp = rbind(Xe[!b.rows,],Xep)
  # plot raw & relative incidence
  plot.ep(Xepp,'raw','HIV Incidence (per person-year)')
    fig.save(uid,N$sam,'foi.ep.incidence.raw',w=9,h=6.5)
  plot.ep(Xepp,'1-2/2','Relative Difference in HIV Incidence (X - NP / NP)',leg='none')
    fig.save(uid,N$sam,'foi.ep.incidence.rel',w=9,h=6)
}

main.wiw = function(){
  X = rbind(
    cbind(clean.wiw.data(read.csvs('foi-ep','wiw','foi','all')),par='Equal Parameters'),
    cbind(clean.wiw.data(read.csvs('fit',   'wiw','foi','all')),par='Recalibrated Parameters'))
  g = do.margin(X,'part',type='rel',strat=c('case.lab','par')) +
    facet_grid('par~case.lab') +
    labs(y='Yearly Infections (%)')
  fig.save(uid,N$sam,'foi.wiw.part',w=11,h=5)
}

main.tpaf = function(){
  X = read.csvs('tpaf','expo','foi')
  X$t.hor = X$t - X$tpaf.t0
  # plot groups of transmission pathways
  pops = list(part=c('msp','cas','swx'),popf=c('aqf','fswf','clif'))
  for (pop in names(pops)){
    g = plot.tpaf.box(X,tpaf.pop=pops[[pop]],tpaf.t0=c(1990,2000,2010)) +
      facet_grid('tpaf.t0 ~ tpaf.pop',scales='free_y')
    g = plot.clean.foi(g)
    fig.save(uid,N$sam,'tpaf.foi',pop,w=8,h=7)
  }
}

# main.ep()
# main.wiw()
# main.tpaf()
