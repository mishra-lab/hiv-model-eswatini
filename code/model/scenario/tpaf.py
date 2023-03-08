from utils import log,fio,dict_list_update
from model import system,params,out
from model.scenario import N,tvec,fname

case = 'base'

def run():
  log(0,'scenario.tpaf.run: '+case)
  P1s = fio.load(fname('npy','fit','Ps',case=case))
  R1s = system.run_n(P1s,t=tvec['main'])
  ekwds = dict(R1s=R1s,tvec=tvec['main'],t=tvec['plot'],vsop='1-2/1',
    snames=['all','w','m','aq','fsw','cli'])
  tpafs = dict(
    msp = dict(p=0),
    cas = dict(p=1),
    swx = dict(p=(2,3)),
    aqf = dict(fs=(0,1),fi=(0,1)),
    aqt = dict(ts=(0,1),ti=(0,1)),
    fswf = dict(fs=0,fi=(2,3)),
    fswt = dict(ts=0,ti=(2,3)),
    clif = dict(fs=1,fi=(2,3)),
    clit = dict(ts=1,ti=(2,3)),
  )
  t0s = [1990,1995,2000,2005,2010,2015,2020,2025]
  E = out.expo([],[],[],[],[],ecols={'tpaf.pop':None,'tpaf.t0':None})
  for name,spec in tpafs.items():
    for t0 in t0s:
      log(1,'tpaf: '+name+'_'+str(t0))
      P2s = dict_list_update(P1s,mix_mask_tpaf=params.get_mix_mask(**spec),t0_tpaf=t0)
      R2s = system.run_n(P2s,t=tvec['main'])
      Ei = out.expo(['cuminfect'],**ekwds,R2s=R2s,ecols={'tpaf.pop':name,'tpaf.t0':str(t0)})
      E = {col:E[col]+Ei[col] for col in E}
  fio.save_csv(fname('csv','tpaf','expo',case=case),E)
