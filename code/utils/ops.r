library('parallel')
options(width=200,scipen=99)
len = length
pout = function(...){ cat(paste0(...,'\n')) } # cleaner print
root.path = function(...,create=FALSE){
  # like fio.rootpath
  root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]
  path = file.path(root,...)
  if (create & !dir.exists(dirname(path))){ dir.create(dirname(path),recursive=TRUE) }
  return(path)
}
p.stars = function(p,s=5){
  paste0(rep('*',min(s,max(0,-ceiling(log10(p))))),collapse='')
}
rename.cols = function(X,...){
  # e.g. rename.cols(X,a='b') renames existing column 'a' to new column 'b'
  map = list(...)
  for (name in names(map)){
    X[[map[[name]]]] = X[[name]]
    X[[name]] = NULL
  }
  return(X)
}
filter.cols = function(X,...){
  # e.g. filter.cols(X,a='b') selects rows in X where X$a %in% c('b')
  filter = list(...)
  bool = !logical(nrow(X))
  for (col in names(filter)){ bool = bool & X[[col]] %in% filter[[col]] }
  X[bool,]
}
rescale = function(x){
  (x - min(x)) / (max(x) - min(x))
}
iqr = function(x){
  unname(diff(quantile(x,c(.25,.75))))
}
par.lapply = function(...,cores=7,.par=TRUE){
  # simple wrapper for parallel::mclapply with some default arguments
  if (.par){
    mclapply(...,mc.cores=cores,mc.set.seed=FALSE)
  } else {
    lapply(...)
  }
}
rbind.lapply = function(...,.par=TRUE){
  # par.lapply & rbind the result
  do.call(rbind,par.lapply(...,.par=.par))
}
