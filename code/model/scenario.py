import numpy as np
from utils import rootpath,genpath,flatten,minimize,itslice,squarish
from utils import stats,fio,deco,parallel
from model import slicers,params,system,target,out,plot,fit

plotsize = 3 # inches
uid = fio.datestamp()
low = (.60,.80,.80)

def fname(case,key,N,N0):
  if case=='npy':
    path,ext = ['data','.npy',uid],''
  if case=='csv':
    path,ext = ['data','mid',uid],'.csv'
  if case=='fig':
    path,ext = ['out','fig',uid],'.pdf'
  return genpath(rootpath(*path,'{}_N={}-{}{}'.format(key,N0,N0+N-1,ext)))

def get_sample(t,T,N,N0=0,top=.10):
  Ps_sam = params.get_n_all(N,seeds=range(N0,N0+N))
  Rs_sam = system.run_n(Ps_sam,t,T)
  Rs_sam = system.drop_fails(Rs_sam)[0]
  fit.plot_all(t,Rs_sam,T,fname=fname('fig','fit_sam',N,N0))
  Rs = target.top_q_ll(Rs_sam,top=top)
  return Rs

def Rqx_by_group(Rqx=None,**kwds):
  # create/update a relative-rate using pop=value kwds, e.g. FSW=1.5
  if Rqx is None: Rqx = np.ones((2,4,1))
  for sname,sRqx in kwds.items():
    for sk in flatten(slicers[sname].pop.get('s',(0,1))):
      for ik in flatten(slicers[sname].pop.get('i',(0,1,2,3))):
        Rqx[sk,ik,0] = sRqx
  return Rqx

def fit_cascade_jfun(pfits,P,PD,T,t):
  # map pfits (quantiles) -> values from PD & update P with the new values
  for pfit,key in zip(pfits,PD.keys()):
    rate,pop = key.split(':')
    P.update({rate+'_scen':Rqx_by_group(Rqx=P.get(rate),**{pop:PD[key].ppf(pfit)})})
  # run the model
  R = system.run(P,t,T,RPts=[])
  return -R['ll']

def fit_cascade(P,PD,T,t,ftol=.1):
  # TODO: re-implement specified P0
  # fit the model to "T", with base parameters P, and re-sampled values in PD
  # PD has keys like "rate:pop"; we re-sample based on percentiles (pfit)
  kwds = dict(method='SLSQP',options=dict(ftol=ftol))
  # we fit relative rates stepwise (dx, +tx, +ux) as it is more efficient
  PDz = {}; P0z = []; Tz = []; # init step-wise: distrs, x0, targets
  for pz,oz in zip(['Rdx','Rtx','Rux'],['diagnosed','treated','vls']):
    PDzi = {k:v for k,v in PD.items() if (pz in k)}
    PDz.update(PDzi)
    P0z += [.1 for k in PDzi]
    Tz  += [Ti for Ti in T if (oz in Ti.name)]
    jfun = lambda pfits: fit_cascade_jfun(pfits,P,PDz,Tz,t)
    M = minimize(jfun,P0z,bounds=[(0,1) for k in PDz],**kwds)
    P0z[:] = M.x
    if not M.success:
      return False
  # DEBUG: non-stepwise
  # jfun = lambda pfits: fit_cascade_jfun(pfits,P,PD,T,t)
  # M = minimize(jfun,[.1 for k in PD],bounds=[(0,1) for k in PD],**kwds)
  return P # TODO: is this robust updated P?

def fit_cascade_n(Ps,PD,T,t,**kwds):
  fun = lambda P: fit_cascade(P,PD=PD,T=T,t=t,**kwds)
  return parallel.ppool(len(Ps)).map(fun,Ps)

def get_sens_data(R,t,case):
  # extract outputs of interest & selected model params for sensitivity analysis
  to  = list(range(2000,2050+1,5))
  ito = itslice(to,t)
  Pkeys = ['seed',
    'EHY_acute','Rbeta_vi_rec', # HIV
    'RPA_condom_a:v','PA_circum_2050', # sex
    'RP_gud_2050','Rbeta_gud_sus_w','Rbeta_gud_inf', # GUD
    'PX_fsw','dur_fsw_l','A_swq_cli','dur_cli', # FSW & Cli
  ]
  D = dict(case=case,
    **{'cuminf_'+str(toi):oi for toi,oi in zip(to,out.cuminfect(R,**slicers['ALL'].pop,aggr=True,tvec=t)[ito])},
    **{'incid_' +str(toi):oi for toi,oi in zip(to,out.incidence(R,**slicers['ALL'].pop,aggr=True,tvec=t,t=to))},
    **{key:R['P'][key] for key in Pkeys},
    **{
      'prev_all_2020': out.prevalence(R,**slicers['ALL'].pop,aggr=True,tvec=t,t=2020)[0],
      'prev_ratio_fsw:wq': out.vs_pop('prevalence',R,slicers['FSW'].pop,slicers['WQ'].pop,tvec=t,t=2020)[0],
      'prev_ratio_cli:mq': out.vs_pop('prevalence',R,slicers['Cli'].pop,slicers['MQ'].pop,tvec=t,t=2020)[0],
    }
  )
  return D

def get_sens_data_n(Rs,t,case):
  return [get_sens_data(R,t,case) for R in Rs]

def plot_diff(t,R1s,R2s,oname,snames,vsop,intervals=.95,
    ylim=None,xlim=None,ylab=None,**kwds):
  assert len(R1s)==len(R2s), 'R1s and R2s must have the same length'
  R1s,R2s = system.drop_fails(R1s,R2s)
  fh,ah = plot.subplots(1,len(snames))
  if ylab is None: ylab = out.labels.get(oname,oname)
  for i,sname in enumerate(snames): # subplots
    plot.plt.sca(ah[0,i])
    plot.labels(title=slicers[sname].label,
      x='Time (years)',y=ylab if i==0 else None)
    for interval in flatten(intervals):
      plot.plot_SvR(oname,t,R1s,R2s,sname,vsop=vsop,tvec=t,box=5,**kwds)
    plot.lims(y=ylim,x=xlim)
    plot.plt.grid(axis='y',lw=.5,zorder=0,color=(.8,.8,.8))
  fh.set_size_inches((plotsize*len(snames),plotsize))
  fh.tight_layout()
  return fh,ah

def get_refit_case(case):
  eps = 1e-9
  PDD = {
    '0-1':  stats.uniform(l=0,h=1),
    '1-1':  stats.uniform(l=1-eps,h=1+eps),
    '1-50': stats.uniform(l=1,h=50),
  }
  if case=='LoLo':
    T2 = target.make_targets_2020(low,[100,100,100],s=(0,1),i=(0,1,2,3))\
       + target.make_targets_2020(low,[100,100,100],s=0,i=(2,3))
    PD = {'Rdx:ALL':PDD['0-1'],'Rtx:ALL':PDD['0-1'],'Rux:ALL':PDD['1-50'],
          'Rdx:FSW':PDD['0-1'],'Rtx:FSW':PDD['0-1'] }
  if case=='LoHi':
    T2 = target.make_targets_2020(low,[100,100,100],s=(0,1),i=(0,1,2,3))
    PD = {'Rdx:ALL':PDD['0-1'],'Rtx:ALL':PDD['0-1'],'Rux:ALL':PDD['1-50'],
          'Rdx:FSW':PDD['1-1'],'Rtx:FSW':PDD['1-1'] }
  if case=='HiLo':
    T2 = target.make_targets_2020(low,[100,100,100],s=0,i=(2,3))
    PD = {'Rdx:ALL':PDD['1-1'],'Rtx:ALL':PDD['1-1'],'Rux:ALL':PDD['1-1'],
          'Rdx:FSW':PDD['0-1'],'Rtx:FSW':PDD['0-1'],'Rux:FSW':PDD['1-50'] }
  return T2,PD

def main(N,N0=0,sample=True,cf=True,refit=True,top=.10,sens=True,plotfit=True,sankey=False):
  t  = system.f_t(t1=2025)
  t2 = system.f_t(t1=2050,dt=.05)
  T1 = target.get_all_esw()
  # base case (fitted)
  if sample:
    R1s = get_sample(t,T1,N,N0=N0,top=top)
    P1s = fio.save(fname('npy','P1s',N,N0),[R['P'] for R in R1s])
  else:
    P1s = fio.load(fname('npy','P1s',N,N0))
  # re-run for t2
  R1s = system.run_n(P1s,t2,T1)
  if sens: fio.save_csv(fname('csv','sens_base',N,N0),get_sens_data_n(R1s,t2,'base'))
  if sankey: fio.save_csv(fname('csv','sankey_base',N,N0),out.get_sankey(R1s,t2,range(1985,2050)))
  if plotfit: fit.plot_all(t2,R1s,T1,fname=fname('fig','fit_base',N,N0))
  if not cf: return
  # counterfactuals
  for case in ['LoLo','LoHi','HiLo']:
    T2,PD = get_refit_case(case)
    casestr = case+'_'+''.join(str(int(z*100)) for z in low)
    if refit:
      P2s = fit_cascade_n(P1s,PD,T2,t,ftol=.1)
      fio.save(fname('npy','P2s_'+casestr,N,N0),P2s)
    else:
      P2s = fio.load(fname('npy','P2s_'+casestr,N,N0))
    # re-run for t2
    R2s = system.run_n(P2s,t2,T2)
    if sens: fio.save_csv(fname('csv','sens_'+casestr,N,N0),get_sens_data_n(R2s,t2,case))
    if sankey: fio.save_csv(fname('csv','sankey_'+casestr,N,N0),out.get_sankey(R1s,t2,range(1985,2050)))
    if plotfit: fit.plot_all(t2,R2s,T2,fname=fname('fig','fit_'+casestr,N,N0))
    # plot cumulative infections averted / extra
    plot_diff(t2,R2s,R1s,'cuminfect',['ALL','FSW','Cli'],intervals=[.5,.95],ylim=(0,.8),
      xlim=(1998,2042),vsop='1-2/1',ylab='Cumulative Infections Averted')
    plot.save(fname('fig','ci_avert_'+casestr,N,N0))
    plot_diff(t2,R2s,R1s,'cuminfect',['ALL','FSW','Cli'],intervals=[.5,.95],ylim=(0,.8),
      xlim=(1998,2042),vsop='1-2/2',ylab='Cumulative Additional Infections')
    plot.save(fname('fig','ci_extra_'+casestr,N,N0))

