import numpy as np
from copy import deepcopy
from utils import _,log,tarray,fio
from model import system,out,fit
from model.scenario import fname

tf = 2050
tlines = { # timelines
  'opt': {'recov':.90,'t5':np.append(np.cumsum([2025,0.25,0.25,0.50]),tf)},
  'pes': {'recov':.60,'t5':np.append(np.cumsum([2025,0.25,1.00,5.00]),tf)},
  'non': {'recov':.00,'t5':np.append(np.cumsum([2025,0.25,1.00,5.00]),tf)},
}
servs = { # service disruption scenarios
  'px':   ['cdm'],
  'tx':   ['dx','unvx','revx'],
  'ax':   ['cdm','dx','unvx','revx'],
}
tvec = { # time vectors
  'main': system.get_t(tf=tf),
  'plot': system.get_t(tf=tf,dt=0.2),
}
ekwds = dict( # keywords for out.expo
  tvec   = tvec['main'],
  t      = tvec['plot'],
  skeys  = ['all','w','m','aq','fsw','cli'],
  onames = ['incidence','prevalence','cuminfect','cumdeath',
    'diagnosed','treated_c','treated_u','vls_c','vls_u'],
  mode = 'id')
# custom fit plots
fit.specsets.update(future=['prev','inc','condom','dx_rate','tx_rate',
  'diag','treat_c','treat_u','vls_c','vls_u'])
for k in fit.specsets['future']: fit.specs[k].update(xlim=(1980,2035))

def edit_tarray(Xo,te,xe):
  b = Xo.ti < te[0] # bool: which ti,xi to keep
  if xe.shape[-1] == len(te)-1: # prepend Xo(te[0]) to Xe
    xe = np.concatenate((Xo(te[0])[...,_],xe),axis=-1)
  return tarray.tarray( # new tarray from (ti,te) and (xi,Xe)
    np.concatenate((Xo.ti[b],te),axis=0),
    np.concatenate((Xo.xi[...,b],xe),axis=-1))

def apply_RR4(x4,RR,recov):
  RR = np.reshape(RR,x4[...,0].shape)
  x4[...,-4] *= RR
  x4[...,-3] *= RR
  x4[...,-2] *= (RR + (1-RR)*recov)
  x4[...,-1] *= (RR + (1-RR)*recov)
  return x4

def get_scen(Ps,adjs,t5,recov=1):
  for P in Ps:
    if 'cdm' in adjs:
      x4 = P['PF_condom_t']([t5[0]]*4)
      x4 = apply_RR4(x4,[.75,.50,.25,.50],recov)
      P['PF_condom_t'] = edit_tarray(P['PF_condom_t'],t5,x4)
    # TODO: PF_circum_t ?
    if 'dx' in adjs:
      x4 = P['dx_sit']([t5[0]]*4)
      x4 = apply_RR4(x4,[[.60,.60,.40,.40],[.40,.40,.40,.40]],recov)
      P['dx_sit'] = edit_tarray(P['dx_sit'],t5,x4)
    if 'unvx' in adjs:
      x4 = P['unvx_sit']([t5[0]]*4)
      x4 = apply_RR4(x4,[[1.5,1.5,3.0,3.0],[2.0,2.0,2.0,2.0]],recov)
      P['unvx_sit'] = edit_tarray(P['unvx_sit'],t5,x4)
    if 'revx' in adjs:
      # special case, since revx_t did not support groups (si)
      b = P['revx_t'].ti < t5[0]
      ti = [*P['revx_t'].ti[b],*t5]
      xi = np.concatenate((
        P['revx_t'].xi[b],      # original xi
        P['revx_t']([t5[0]]*5), # repeat (t5[0]) like edit_tarray
      ),axis=-1) * np.ones((2,4,1,1,1))
      xi = apply_RR4(xi,[[0.75,0.75,0.50,0.50],[.50,.50,.50,.50]],recov)
      P['revx_t'] = tarray.tarray(ti,xi)
  return Ps

def run(base=True):
  log(0,'future.run')
  P0s = fio.load_npy(fname('npy','fit','Ps',case='base'))
  if base:
    R0s = system.run_n(P0s,t=tvec['main'])
    fio.save_csv(fname('csv','future','expo',case='base'),
      out.expo(R0s,**ekwds,ecols=dict(serv='base',tline='base')))
    fit.plot_sets(tvec['main'],R0s,tfname=fname('fig','future','{}',case='base'),sets='future')
  for serv,adjs in servs.items():
    for tline,kwds in tlines.items():
      case = serv+tline[0]
      Ps = get_scen(deepcopy(P0s),adjs,**kwds)
      Rs = system.run_n(Ps,t=tvec['main'])
      fio.save_csv(fname('csv','future','expo',case=case),
        out.expo(Rs,**ekwds,ecols=dict(serv=serv,tline=tline)))
      fit.plot_sets(tvec['main'],Rs,tfname=fname('fig','future','{}',case=case),sets='future')

if __name__ == '__main__':
  run()
