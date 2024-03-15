# function to compute TPAFs by running model with vs without masked transmission

from utils import log,dict_list_update
from model import system,params,out

ekwds = dict( # default kwds for out.expo
  skeys = ['all','w','m','aq','asw','fsw','cli'],
  onames = ['cuminfect'],
  vsop = '1-2/1',
  mode = 'q')

def run(Ps,tvec,t,paths,t0s):
  # runs the model for Ps & tvec, then again for all combinations of paths (dict) & t0s (list)
  # and computes the corresponding TPAFs for years t in the list-of-lists E from out.expo
  # e.g. paths = dict(msp=dict(p=0,...), t0 = [2000], t = [2010] computes
  # the 10-year TPAF starting from 2000 named 'msp' defined as in get_mix_mask(p=0)
  ekwds.update(tvec=tvec,t=t)
  R1s = system.run_n(Ps,t=tvec)
  E = out.expo(R1s,[],[],[],[],ecols={'tpaf.path':None,'tpaf.t0':None},mode=ekwds['mode'])
  for name,inds in paths.items():
    for t0 in t0s:
      log(1,'tpaf: '+name+' @ '+str(t0))
      P2s = dict_list_update(Ps,mix_mask_tpaf=params.get_mix_mask(**inds),t0_tpaf=t0)
      R2s = system.run_n(P2s,t=tvec)
      Ei = out.expo(R1s=R1s,R2s=R2s,**ekwds,t0=t0,ecols={'tpaf.path':name,'tpaf.t0':str(t0)})
      E = {col:E[col]+Ei[col] for col in E}
  return E
