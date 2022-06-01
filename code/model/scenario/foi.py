import numpy as np
from copy import copy
from utils import log,fio,itslice,dict_list_update
from model import target,system,fit,out,slicers,params
from model.scenario import tvec,fname
import model.scenario

model.scenario.uid = '2022-06-01'
model.scenario.N['sam'] = 100000

cases = [
  'base',
  'bpd',
  'bpy',
  'bmy',
]

def get_keyout_data(R,t,case):
  pops = ['all','w','m','aq','wq','mq','fsw','cli']
  to = tvec['outs'].tolist()
  Pkeys = ['seed','foi_mode'] + list(params.def_sample_distrs().keys()) + ['PX_fsw','PX_cli','EHY_acute']
  return dict(case=case,ll=R.get('ll',np.nan),
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
  Ps = fio.load(fname('npy','fit','Ps',case=case))
  Rs = system.run_n(Ps,t=tvec['main'],T=T)
  fit.plot_cal(tvec['main'],Rs,T,fname=fname('fig','foi','cal',case=case))
  fio.save_csv(fname('csv','foi-fit','keyout',case=case),
    [get_keyout_data(R,tvec['main'],case) for R in Rs])
  fio.save_csv(fname('csv','foi-fit','infs',case=case),
    out.get_infections(Rs,tvec['main'],tvec['plot']))

def run_ep(top=1.):
  log(0,'scenario.foi.run_ep: '+', '.join(cases))
  onames = ['incidence','prevalence']
  ekwds = dict(tvec=tvec['main'],t=tvec['plot'],snames=['all','w','m','aq','cli','fsw'])
  Ps = target.top_ll(fio.load(fname('npy','fit','Ps',case='base')),top)
  for case in cases:
    log(1,case)
    R1s = system.run_n(dict_list_update(Ps,foi_mode=case),t=tvec['main'])
    fio.save_csv(fname('csv','foi-ep','keyout',case=case),
      [get_keyout_data(R,tvec['main'],case) for R in R1s])
    fio.save_csv(fname('csv','foi-ep','infs',case=case),
      out.get_infections(R1s,tvec['main'],tvec['plot']))
    E = out.expo(onames,R1s,**ekwds)
    if case == 'base':
      R2s = copy(R1s)
    else:
      E_vs = out.expo(onames,R1s,R2s=R2s,vsop='1-2',**ekwds)
      E = {k:E[k]+E_vs[k] for k in E}
    fio.save_csv(fname('csv','foi-ep','expo',case=case),E)

def run_tpaf(case,top=1.):
  log(0,'scenario.foi.run_tpaf: '+case)
  P1s = target.top_ll(fio.load(fname('npy','fit','Ps',case=case)),top)
  R1s = system.run_n(P1s,t=tvec['main'])
  ekwds = dict(R1s=R1s,tvec=tvec['main'],t=tvec['plot'],snames=['all','w','m','aq','cli','fsw'],vsop='1-2/1')
  tpafs = dict(
    msp   = dict(p=0),
    cas   = dict(p=1),
    swq   = dict(p=(2,3)),
    aqfr  = dict(sfr=(0,1),ifr=(0,1)),
    fswfr = dict(sfr=0,ifr=(2,3)),
    clifr = dict(sfr=1,ifr=(2,3)),
    aqto  = dict(sto=(0,1),ito=(0,1)),
    fswto = dict(sto=0,ito=(2,3)),
    clito = dict(sto=1,ito=(2,3)),
  )
  t0s = [1990,1995,2000,2005,2010,2015,2020,2025,2030]
  E = out.expo([],[],[],[],[],ecols=dict(tpaf=None))
  for tpaf,spec in tpafs.items():
    for t0 in t0s:
      log(1,tpaf+'_'+str(t0))
      P2s = dict_list_update(P1s,mix_mask_tpaf=params.get_mix_mask(**spec),t0_tpaf=t0)
      R2s = system.run_n(P2s,t=tvec['main'])
      E1 = out.expo(['cuminfect'],**ekwds,R2s=R2s,ecols={'tpaf.pop':tpaf,'tpaf.t0':str(t0)})
      E = {k:E[k]+E1[k] for k in E}
  fio.save_csv(fname('csv','foi-tpaf','expo',case=case),E)

if __name__ == '__main__':
  # run_fit()
  # run_ep()
  # run_tpaf()
  pass
