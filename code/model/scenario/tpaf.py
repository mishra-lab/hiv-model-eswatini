from utils import log,fio,dict_list_update
from model import system,params,out
from model.scenario import tvec,fname

tkp = dict(tvec=tvec['main'],t=tvec['plot'])
ekwds = dict(
  snames = ['all','w','m','aq','asw','fsw','cli'],
  onames = ['cuminfect'],
  vsop = '1-2/1',
  mode = 'q')
tpafs = dict(
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
)
t0s = [1990,2000,2010,2020]

def run(case='base'):
  log(0,'scenario.tpaf.run: '+case)
  P1s = fio.load(fname('npy','fit','Ps',case=case))
  R1s = system.run_n(P1s,t=tvec['main'])
  E = out.expo([],[],[],[],[],ecols={'tpaf.pop':None,'tpaf.t0':None})
  for name,spec in tpafs.items():
    for t0 in t0s:
      log(1,'tpaf: '+name+'_'+str(t0))
      P2s = dict_list_update(P1s,mix_mask_tpaf=params.get_mix_mask(**spec),t0_tpaf=t0)
      R2s = system.run_n(P2s,t=tvec['main'])
      Ei = out.expo(R1s=R1s,R2s=R2s,**ekwds,**tkp,t0=t0,ecols={'tpaf.pop':name,'tpaf.t0':str(t0)})
      E = {col:E[col]+Ei[col] for col in E}
  # fio.save_csv(fname('csv','tpaf','expo',case=case),E)
