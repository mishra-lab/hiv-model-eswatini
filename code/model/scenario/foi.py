from copy import copy
from utils import log,fio,dict_list_update
from model import system,out
from model.scenario import calibrate,tpaf,tvec,fname
import model.scenario

model.scenario.uid = '2023-03-03'
# model.scenario.N['sam'] = 1000 # DEBUG

cases = ['base','rd','ry','pd','py']

kwds = fio.argvkwds(case='base',scinet=False)
if kwds.pop('scinet'): model.scenario.scinet()

def eprun():
  log(0,'foi.eprun: '+', '.join(cases))
  onames = ['incidence','prevalence']
  ekwds = dict(tvec=tvec['main'],t=tvec['plot'],snames=['all','w','m','aq','fsw','cli'])
  Ps = fio.load(fname('npy','fit','Ps',case='base'))
  for case in cases:
    log(1,case)
    R1s = system.run_n(dict_list_update(Ps,foi_mode=case),t=tvec['main'])
    E = out.expo(onames,R1s,**ekwds)
    fio.save_csv(fname('csv','foi-ep','wiw',case=case),out.wiw(R1s,tvec['main'],tvec['plot']))
    if case == 'base':
      R2s = copy(R1s)
    else:
      EAD = out.expo(onames,R1s,R2s=R2s,vsop='1-2',  **ekwds)
      ERD = out.expo(onames,R1s,R2s=R2s,vsop='1-2/2',**ekwds)
      E = {col:E[col]+EAD[col]+ERD[col] for col in E}
    fio.save_csv(fname('csv','foi-ep','expo',case=case),E)

if __name__ == '__main__':
  # calibrate.run(**kwds)
  # calibrate.merge(**kwds)
  # calibrate.rerun(**kwds)
  # tpaf.run(**kwds)
  # eprun()
  pass
