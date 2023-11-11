from copy import copy
from utils import log,fio,dict_list_update
from model import system,out
from model.scenario import imis,tpaf,akwds,tvec,fname

cases = ['base','foi-rd','foi-ry','foi-py']

tkp = dict(tvec=tvec['main'],t=tvec['plot'])
ekwds = dict(
  snames = ['all','w','m','aq','fsw','cli'],
  onames = ['incidence','prevalence'],
  mode = 'q')

def run_ep():
  log(0,'foi.run_ep: '+', '.join(cases))
  Ps = fio.load(fname('npy','fit','Ps',case='base'))
  for case in cases:
    log(1,case)
    R1s = system.run_n(dict_list_update(Ps,foi_mode=case),t=tvec['main'])
    E = out.expo(R1s,**tkp,**ekwds)
    fio.save_csv(fname('csv','foi-ep','wiw',case=case),out.wiw(R1s,**tkp))
    if case == 'base':
      R2s = copy(R1s)
    else:
      EAD = out.expo(R1s,R2s=R2s,vsop='1-2',  **tkp,**ekwds)
      ERD = out.expo(R1s,R2s=R2s,vsop='1-2/2',**tkp,**ekwds)
      E = {col:E[col]+EAD[col]+ERD[col] for col in E}
    fio.save_csv(fname('csv','foi-ep','expo',case=case),E)

if __name__ == '__main__':
  # imis.run(**akwds,foi_mode=akwds['case'][4:])
  # akwds.pop('b')
  # imis.sample_post(**akwds)
  # imis.rerun(**akwds)
  # tpaf.run(**akwds)
  # run_ep()
  pass
