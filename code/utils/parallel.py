from pathos.multiprocessing import ProcessingPool

global cpus

def ppool(N,*args,**kwds):
  global cpus
  return ProcessingPool(min(N,cpus),*args,**kwds)