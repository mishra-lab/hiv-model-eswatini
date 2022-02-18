options(width=200)
root.path = function(...,create=FALSE){
  root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]
  path = file.path(root,...)
  if (create & !dir.exists(dirname(path))){ dir.create(dirname(path),recursive=TRUE) }
  return(path)
}
quiet = function(code){
  sink('/dev/null')
  result = code
  sink()
  return(result)
}
str.r.pad = function(s,n){
  return(sprintf(fmt=paste0('%-',n,'s'),s))
}
str.l.pad = function(s,n){
  return(sprintf(fmt=paste0('%+',n,'s'),s))
}
print.v = function(v,w=8,n=3){
  return(sprintf(fmt=paste0('%',w,'.',n,'f'),v))
}
print.m.ci = function(m,ci,w=8,n=3){
  return(paste0('',print.v(m,w=w,n=n),
    ' (',print.v(ci[1],w=w,n=n),',',print.v(ci[2],w=w,n=n),')'))
}
print.p = function(p,w=7,n=5,s=5){
  stars = rep('*',min(s,max(0,-ceiling(log10(p)))))
  return(paste0(print.v(p,w=w,n=n),' ',str.r.pad(paste0(stars,collapse=''),s)))
}
squarish = function(n){
  i = ceiling(sqrt(n))
  j = ceiling(n/i)
  return(c(i,j))
}
iterms = function(x,n,lower=TRUE){
  return(paste(lapply(ifelse(lower,1,n):n,function(ni){
    paste(apply(combn(x,ni),2,paste,collapse=':'),collapse=' + ')
  }),collapse=' + '))
}
block.sample = function(X,id.var='id',n=NULL,replace=TRUE){
  idu = unique(X[[id.var]])
  if (missing(n)){ n = length(idu) }
  ids = sample(idu,n,replace=replace)
  return(do.call(rbind,lapply(ids,function(id){
    return(X[X[[id.var]]==id,])
  })))
}