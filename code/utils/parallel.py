from pathos.multiprocessing import ProcessingPool

cpus = 7

def ppool(N,*args,**kwds):
  global cpus
  return ProcessingPool(min(N,cpus),*args,**kwds)