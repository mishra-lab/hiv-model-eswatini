import numpy as np
from utils import dict_split
from model import slicers,out

def check(R,sname1,sname2,onames=None,tol=1e-9,**kwds):
  if onames is None:
    onames = ['prevalence','incidence','diagnosed','treated_c','vls_c']
  fkwds = dict_split(kwds,['tvec'])
  S1 = slicers[sname1]
  S2 = slicers[sname2]
  for oname in onames:
    ofun = out.by_name(oname)
    o1 = ofun(R,**S1.pop,**fkwds,)
    o2 = ofun(R,**S2.pop,**fkwds,)
    ok = np.allclose(o1,o2,rtol=1e-9,atol=1e-12)
    print('collapse: {} vs {} @ {}: {}'.format(
      sname1.rjust(6),sname2.ljust(6),oname.rjust(12),'OK' if ok else 'FAIL (!)'))
