library('parallel')
options(width=200,scipen=99)
len = length
root.path = function(...,create=FALSE){
  root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]
  path = file.path(root,...)
  if (create & !dir.exists(dirname(path))){ dir.create(dirname(path),recursive=TRUE) }
  return(path)
}
read.big.csv = function(fname,fresh=FALSE,...){
  # for big "path/file.csv", write "path/.tmp/file.csv.rds" on first run & use this thereafter
  fname.rds = file.path(dirname(fname),'.tmp',paste0(basename(fname),'.rds'))
  if (fresh){
    X = read.csv(fname,...)
    if (!dir.exists(dirname(fname.rds))) { dir.create(dirname(fname.rds)) }
    saveRDS(X,fname.rds)
    return(X)
  } else {
    return(readRDS(fname.rds))
  }
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
p.stars = function(p,s=5){
  paste0(rep('*',min(s,max(0,-ceiling(log10(p))))),collapse='')
}
p.print = function(p,w=7,n=5,s=5){
  paste0(print.v(p,w=w,n=n),' ',str.r.pad(p.stars(p,s),s))
}
squarish = function(n){
  i = ceiling(sqrt(n))
  j = ceiling(n/i)
  return(c(i,j))
}
rename.cols = function(X,...){
  # e.g. rename.cols(X,a=b) renames existing column 'a' to new column 'b'
  map = list(...)
  for (name in names(map)){
    X[[map[[name]]]] = X[[name]]
    X[[name]] = NULL
  }
  return(X)
}
filter.cols = function(X,...){
  filter = list(...)
  bool = !logical(nrow(X))
  for (col in names(filter)){ bool = bool & X[[col]] %in% filter[[col]] }
  X[bool,]
}
iqr = function(x){
  unname(diff(quantile(x,c(.25,.75))))
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
.n.cores = 7
par.lapply = function(...,cores=.n.cores,.par=TRUE){
  # simple wrapper for parallel::mclapply with some default arguments
  if (.par){
    mclapply(...,mc.cores=cores,mc.set.seed=FALSE)
  } else {
    lapply(...)
  }
}
