source('post/config.r')
source('post/post.r')
source('post/wiw.r')
source('post/tpaf.r')

plot.clean.foi = function(g,leg='top',lab='FOI Approach'){
  g = g + labs(color=lab,fill=lab,linetype=lab) +
    scale_fill_manual(values=set.cols$foi) +
    scale_color_manual(values=set.cols$foi) +
    scale_linetype_manual(values=set.lts$foi) +
    theme(legend.position=leg)
}

main.post = function(){
  X = read.csvs('fit','Ps','foi',rdata='load')
  g = plot.post.uni(X,'case.lab',ncol=7,stats=TRUE) # MAN
  g = plot.clean.foi(g)
  fig.save(uid,nid,'post.distr.foi',w=12,h=16)
}

main.ll = function(){
  X = read.csvs('fit','Ps','foi')
  X$case.facet = factor(gsub(' ','\n',X$case.lab),gsub(' ','\n',levels(X$case.lab)))
  g = plot.ll.hist(X,fill=case.lab,ll.min=-500) +
    facet_grid('case.facet')
  g = plot.clean.foi(g,leg='none')
  fig.save(uid,nid,'ll.hist.foi',w=4,h=6)
}

plot.ep = function(Xepp,op,ylab,leg='top',oname='incidence'){
  g = plot.expo.ribbon(filter.cols(Xepp,op=op),oname,grp='case.lab') +
    facet_grid('pop ~ case.lab.x',scales='free_y') +
    labs(x='Year',y=ylab)
  g = plot.clean.foi(g,leg=leg)
}

main.ep = function(tmax=2025){
  # load & clean data
  X = read.csvs('foi-ep','expo','foi')
  Xe = filter.cols(X,pop=c('aq','fsw','cli'),t=seq(1980,tmax))
  # Xe = filter.cols(X,pop=c('aq','fsw'),t=seq(1980,2035)) # slides
  Xe$pop = factor(Xe$pop,levels=names(strat.labs),labels=strat.labs)
  Xe$case.lab.x = Xe$case.lab
  # duplicate base data for each other case
  b.cols = grep('^q.*|case.lab$|^op$',names(Xe))
  b.rows = Xe$case=='base'
  Xep = Xe[!b.rows,]
  Xep[b.cols] = Xe[b.rows,b.cols]
  Xepp = rbind(Xe[!b.rows,],Xep)
  # plot raw & relative incidence
  plot.ep(Xepp,'raw','HIV Incidence (per person-year)')
    fig.save(uid,nid,'foi.ep.incidence.raw',w=9,h=6.5)
    # fig.save(uid,nid,'foi.ep.incidence.slides',w=7,h=4) # slides
  plot.ep(Xepp,'1-2/2','Relative Difference in HIV Incidence (X - EPA / EPA)\n',leg='none')
    fig.save(uid,nid,'foi.ep.incidence.rel',w=9,h=6)
}

main.wiw = function(tmax=2025){
  X = rbind(
    cbind(clean.wiw.data(read.csvs('foi-ep','wiw','foi'),tmax=tmax),par='Equal Parameters'),
    cbind(clean.wiw.data(read.csvs('fit',   'wiw','foi'),tmax=tmax),par='Recalibrated Parameters'))
  g = wiw.margin(X,'ptr',type='rel',strat=c('case.lab','par')) +
    facet_grid('par~case.lab') +
    labs(y='Yearly Infections (%)')
  fig.save(uid,nid,'foi.wiw.ptr',w=11,h=5)
}

main.tpaf = function(){
  X = read.csvs('tpaf','expo','foi')
  X$t.hor = X$t - X$tpaf.t0
  # plot groups of transmission pathways
  paths = c('msp','cas','swx','aqf','fswf','clif')
  g = plot.tpaf(X,tpaf.path=paths,tpaf.t0=seq(1990,2020,10),t.hor=3,t='tpaf.t0') +
    facet_wrap('tpaf.path') + labs(x='Year TPAF starts',y='3-Year TPAF (%)')
  g = plot.clean.foi(g)
  fig.save(uid,nid,'foi.tpaf.3x2',w=8,h=5)
}

# main.post()
# main.ll()
# main.ep()
# main.wiw()
# main.tpaf()
