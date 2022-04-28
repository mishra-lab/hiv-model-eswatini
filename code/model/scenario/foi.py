import numpy as np
from utils import log,fio,itslice
from model import target,system,fit,out,slicers,params
from model.scenario import tvec,fname
import model.scenario

model.scenario.uid = '2022-04-20'
model.scenario.N['cal'] = 100000

cases = dict(
  # bpd  = dict(mode='bpd',plotkwds=dict(color=(0.267,0.005,0.329),ls=':', label='$\\langle1a\\rangle$')),
  bpy  = dict(mode='bpy',plotkwds=dict(color=(0.23 ,0.322,0.546),ls='-.',label='$\\langle1b\\rangle$')),
  bmy  = dict(mode='bmy',plotkwds=dict(color=(0.128,0.567,0.551),ls='--',label='$\\langle2b\\rangle$')),
  # lin  = dict(mode='lin',plotkwds=dict(color=(0.369,0.789,0.383),ls='-', label='$\\langle3\\rangle$')),
  base = dict(mode='fpe',plotkwds=dict(color=(0.993,0.906,0.144),ls='-', label='$\\langle4\\!*\\!\\rangle$')),
)

def get_expo_data(R,t,case):
  pops = ['all','w','m','aq','wq','mq','fsw','cli']
  to = tvec['outs'].tolist()
  Pkeys = ['seed','foi_mode'] + list(params.def_sample_distrs().keys()) + ['PX_fsw','PX_cli','EHY_acute']
  return dict(case=case,ll=R['ll'],
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
  fio.save_csv(fname('csv','fit','expo',case=case,b='all'),
    [get_expo_data(R,tvec['main'],case) for R in Rs])
  fio.save_csv(fname('csv','fit','infs',case=case,b='all'),
    out.get_infections(Rs,tvec['main'],tvec['infs']))

run_fit('base')

