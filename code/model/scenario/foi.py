from copy import copy
from utils import log,fio,dict_list_update
from model import system,out,tpaf
from model.scenario import tvec,fname

cases = ['base','foi-rd','foi-ry','foi-py']

tkp = dict(tvec=tvec['main'],t=tvec['plot'])
ep_ekwds = dict(
  skeys = ['all','w','m','aq','fsw','cli'],
  onames = ['incidence','prevalence'],
  mode = 'q')
tpaf_kwds = dict(
  paths = dict(
    msp = dict(p=0),
    cas = dict(p=1),
    swo = dict(p=2),
    swr = dict(p=3),
    swx = dict(p=(2,3)),
    aqf = dict(fs=(0,1),fi=(0,1)),
    aqt = dict(ts=(0,1),ti=(0,1)),
    aswf = dict(fs=(0,1),fi=(2,3)),
    aswt = dict(ts=(0,1),ti=(2,3)),
    fswf = dict(fs=0,fi=(2,3)),
    fswt = dict(ts=0,ti=(2,3)),
    clif = dict(fs=1,fi=(2,3)),
    clit = dict(ts=1,ti=(2,3)),
  ),
  t0s = [1990,2000,2010,2020],
)

def run_ep():
  log(0,'foi.run_ep: '+', '.join(cases))
  Ps = fio.load_npy(fname('npy','fit','Ps',case='base'))
  for case in cases:
    log(1,case)
    R1s = system.run_n(dict_list_update(Ps,foi_mode=case.replace('foi-','')),t=tvec['main'])
    E = out.expo(R1s,**tkp,**ep_ekwds)
    fio.save_csv(fname('csv','foi-ep','wiw',case=case),out.wiw(R1s,**tkp))
    if case == 'base':
      R2s = copy(R1s)
    else:
      Eda = out.expo(R1s,R2s=R2s,vsop='1-2',  **tkp,**ep_ekwds)
      Edr = out.expo(R1s,R2s=R2s,vsop='1-2/2',**tkp,**ep_ekwds)
      E = {col:E[col]+Eda[col]+Edr[col] for col in E}
    fio.save_csv(fname('csv','foi-ep','expo',case=case),E)

def run_tpaf(case):
  log(0,'foi.run_tpaf: '+case)
  P1s = fio.load_npy(fname('npy','fit','Ps',case=case))
  E = tpaf.run(P1s,**tkp,**tpaf_kwds)
  fio.save_csv(fname('csv','tpaf','expo',case=case),E)
