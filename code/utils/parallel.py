from pathos.multiprocessing import ProcessingPool

cpus = 7

def ppool(n,*args,**kwds):
  global cpus
  return ProcessingPool(min(n,cpus),*args,**kwds)