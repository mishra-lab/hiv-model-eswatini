from utils import fio
from model import system,params,plot,strats

t = system.get_t(t0=2000,tf=2030)
P = params.get_n_all([0])[0]
Rtx_ht = P['Rtx_ht'](t)[0,0,0,:,:]

fh,ah = plot.subplots(1,1)
keys = ('ahi','hiv.1','hiv.2','hiv.3','aids')
for h,key in enumerate(keys):
  S = strats[key]
  plot.line(t,Rtx_ht[h,],color=S.color,label=S.label)
plot.labels(y='Relative ART initiation rate')
ah[0,0].legend()
fh.set_size_inches((6,3))
fh.tight_layout()
plot.save(fio.rootpath('out','fig','toy','art.elig.pdf'))
