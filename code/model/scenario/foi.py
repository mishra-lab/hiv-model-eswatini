import gc
import numpy as np
from copy import copy
from utils import log,fio,itslice,dict_list_update
from model import target,system,fit,out,slicers,params
from model.scenario import tvec,fname
import model.scenario

model.scenario.uid = '2022-04-20'
model.scenario.N['cal'] = 100000

cases = dict(
  base = 'fpe',
  bpd  = 'bpd',
  bpy  = 'bpy',
  bmy  = 'bmy',
)

def get_keyout_data(R,t,case):
  pops = ['all','w','m','aq','wq','mq','fsw','cli']
  to = tvec['outs'].tolist()
  Pkeys = ['seed','foi_mode'] + list(params.def_sample_distrs().keys()) + ['PX_fsw','PX_cli','EHY_acute']
  return dict(case=case,ll=R.pop('ll',np.nan),
    **{key:R['P'][key] for key in Pkeys},
    **{'cuminf_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.cuminfect(R,**slicers[pop].pop,tvec=t)[itslice(to,t)])},
    **{'inc_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.incidence(R,**slicers[pop].pop,tvec=t,t=to))},
    **{'prev_'+pop+'_'+str(toi):oi for pop in pops
        for toi,oi in zip(to,out.prevalence(R,**slicers[pop].pop,tvec=t,t=to))},
    **{step+'_'+pop+'_2020':out.by_name(step)(R,**slicers[pop].pop,tvec=t,t=2020)[0] for pop in pops
        for step in ['diagnosed','treated_c','vls_c','treated_u','vls_u']},
  )

def run_fit(case):
  log(0,'scenario.foi.run_fit: '+case)
  T = target.get_all_esw()
  Ps = fio.load(fname('npy','fit','Ps',case=case,b='all'))
  Rs = system.run_n(Ps,t=tvec['main'],T=T)
  fit.plot_cal(tvec['main'],Rs,T,fname=fname('fig','fit','plot',case=case,b='all'))
  fio.save_csv(fname('csv','fit','keyout',case=case,b='all'),
    [get_keyout_data(R,tvec['main'],case) for R in Rs])
  fio.save_csv(fname('csv','fit','infs',case=case,b='all'),
    out.get_infections(Rs,tvec['main'],tvec['plot']))

def run_ep(top=1.):
  log(0,'scenario.foi.run_ep: '+', '.join(cases))
  onames = ['incidence','prevalence']
  ekwds = dict(tvec=tvec['main'],t=tvec['plot'],snames=['all','w','m','aq','cli','fsw'])
  Ps = target.top_ll(fio.load(fname('npy','fit','Ps',case='base',b='all')),top)
  for case in cases:
    log(1,case)
    R1s = system.run_n(dict_list_update(Ps,foi_mode=cases[case]),t=tvec['main'])
    fio.save_csv(fname('csv','foi-ep','keyout',case=case,b='all'),
      [get_keyout_data(R,tvec['main'],case) for R in R1s])
    fio.save_csv(fname('csv','foi-ep','infs',case=case,b='all'),
      out.get_infections(R1s,tvec['main'],tvec['plot']))
    E = out.expo(onames,R1s,**ekwds)
    if case == 'base':
      R2s = copy(R1s)
    else:
      E_vs = out.expo(onames,R1s,R2s=R2s,vsop='1-2',**ekwds)
      E = {k:E[k]+E_vs[k] for k in E}
    fio.save_csv(fname('csv','foi-ep','expo',case=case,b='all'),E)

if __name__ == '__main__':
  # run_fit()
  # run_ep()
  pass
