import re
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib as mpl
np.set_printoptions(suppress=True)

# comparing how per-contact transmission modifiers are applied
# to compute per-partnership probability of transmission (B):
# - BPH: between-partnership heterogeneity (arithmetic.mean) [Bb]
# - WPH: within-partnership heterogeneity (geometric.mean) [Bw]
# beta - per-contact probability
# R    - relative beta for the factor
# a    - (\alpha) proportion of contacts affected by R
# A    - total contacts per partnership
# maxBratio: finds conditions (R,a,A) to maximize Bw / Bb
# maxplot: plots (a,A) which maximize Bw / Bb given R
# surface: big grid of Bw, Bb, and Bratio for pairs of R,A,a

beta = .0034 # Boily2009

def Bb(R,A,a,beta=.0034):
  return (1-(1-R*beta)**A)*(a) + (1-(1-beta)**A)*(1-a)

def Bw(R,A,a,beta=.0034):
  return 1 - ((1-R*beta)**(A*(a)) * (1-beta)**(A*(1-a)))

def Bratio(R,A,a,beta=.0034):
  return Bw(R,A,a,beta) / Bb(R,A,a,beta)

def maxBratio(R0=.5,Rmin=.001,Rmax=100,beta=.0034,log=True):
  x0 = [R0,100,.5] # initial conditions
  xb = [(Rmin,Rmax),(1,10000),(0,1)] # bounds
  jfun = lambda x: -Bratio(*np.power(10,x),beta) # maximize ratio
  M = minimize(jfun,np.log10(x0),method='L-BFGS-B',bounds=np.log10(xb))
  xopt = np.power(10,M.x)
  if log: print('f = {:.3f}\n  R = {:.2f}\n  A = {:0f}\n  a = {:3f}'.format(-M.fun,*xopt)+
    '\nB_bph = {:.4f}\nB_wph = {:.4f}'.format(Bb(*xopt,beta),Bw(*xopt,beta)))
  return (*xopt,-M.fun)

def labs(x,y,t=None):
  plt.xlabel(x)
  plt.ylabel(y)
  plt.title(t)

def sca(ah):
  plt.sca(ah)
  ah.xaxis.grid(True,color='w',alpha=.1,zorder=99)
  ah.yaxis.grid(True,color='w',alpha=.1,zorder=99)

def contour(*args,**kwds):
  # https://stackoverflow.com/questions/8263769
  h = plt.contourf(*args,**kwds)
  for c in h.collections:
    c.set_edgecolor('face')

def cbar(name,lvls,pos=None,labels=None,**kwds):
  if pos is not None:
    kwds.update(cax=plt.gcf().add_axes(pos))
  cb = plt.colorbar(
    plt.cm.ScalarMappable(
      cmap=plt.get_cmap(name,len(lvls)-1),
      norm=mpl.colors.Normalize(lvls[0],lvls[-1])),
    **kwds)
  if cb.orientation == 'vertical':
    cb.ax.set_yticklabels(labels=np.round(labels,3))
  if cb.orientation == 'horizontal':
    cb.ax.set_xticklabels(labels=np.round(labels,3))

def maxplot():
  fh,ah = plt.subplots(1,1,figsize=(6,4))
  cm = plt.cm.get_cmap('viridis',10); cf = '#cccccc'
  Rv = np.round(np.arange(-2,2.02,.02),3)
  for f in (1,0.5,2):
    x0 = 3*[np.nan]
    for R in Rv:
      if abs(R) > 5e-2: x1 = maxBratio(*3*[10**R],beta=f*beta,log=False) # discontinuity bug
      if f == 1:
        kwds = dict(color=cm(2*(x1[3]-1)),lw=4)
        if R == 0:
          plt.text(x1[1]*3,.5,'$\\frac{1}{2}\\beta$',ha='center',va='top',color=cf)
          plt.text(x1[1]/3,.5,'$2\\beta$',ha='center',va='bottom',color=cf)
        if R in [-2,-1,0,+1,+2]:
          plt.semilogx(x1[1],x1[2],'+',zorder=3,color='red')
          plt.text(x1[1],x1[2],re.sub('.0 $',' ','R = {} '.format(10**R)),va='bottom',ha='right',color='red')
      else:
        kwds = dict(color=cf)
      plt.semilogx((x0[1],x1[1]),(x0[2],x1[2]),**kwds)
      x0 = x1
  labs('A','$\\alpha$')
  plt.ylim((0,1))
  plt.xticks(*2*[[10,100,1000,10000]])
  LR = np.linspace(1,1.5,10+1)
  cbar('viridis',LR-1,orientation='vertical',ticks=LR[::2]-1,labels=LR[::2])
  plt.tight_layout()
  plt.savefig('B.xph.max.pdf')

def surface():
  Ni = 32
  R = np.linspace(-2,1,Ni)
  A = np.linspace( 0,3,Ni)
  a = np.linspace( 0,1,Ni)
  LB = np.linspace(-3,0,12+1)
  LR = np.linspace(1,1.5,10+1)
  Rticks = lambda: plt.yticks(ticks=[-2,-1,0,1],labels=[.01,.1,1,10])
  Aticks = lambda: plt.xticks(ticks=[0,1,2,3],labels=[1,10,100,1000])
  plotBw = lambda M1,M2,MR,MA,Ma: contour(M1,M2,np.log10(Bw(10**MR,10**MA,Ma)),levels=LB,extend='both',cmap='inferno')
  plotBb = lambda M1,M2,MR,MA,Ma: contour(M1,M2,np.log10(Bb(10**MR,10**MA,Ma)),levels=LB,extend='both',cmap='inferno')
  plotRB = lambda M1,M2,MR,MA,Ma: contour(M1,M2,Bratio(10**MR,10**MA,Ma),      levels=LR,extend='both',cmap='viridis')
  fh,ah = plt.subplots(4,3,figsize=(9,13))
  tex = dict(BR = '$B_{WPH}~/~B_{BPH}$',
    Bwph = '$B_{WPH} = 1 - \\prod_f~{(1 - R_f\\,\\beta)}^{\\,A\\,\\alpha_f}$',
    Bbph = '$B_{BPH} = 1 - \\sum_f~\\alpha_f\\,{(1 - R_f\\,\\beta)}^{\\,A}$')
  # R vs a
  Ma,MR = np.meshgrid(a,R); MA = 1.5
  sca(ah[0][0]); plotBw(Ma,MR,MR,MA,Ma); Rticks(); labs('$\\alpha$','$R,\\alpha~(A = 32)$\n\n$R$',tex['Bwph']+'\n')
  sca(ah[0][1]); plotBb(Ma,MR,MR,MA,Ma); Rticks(); labs('$\\alpha$',None,tex['Bbph']+'\n')
  sca(ah[0][2]); plotRB(Ma,MR,MR,MA,Ma); Rticks(); labs('$\\alpha$',None,tex['BR']+'\n')
  # R vs A
  MA,MR = np.meshgrid(A,R); Ma = 0.5
  sca(ah[1][0]); plotBw(MA,MR,MR,MA,Ma); Rticks(); Aticks(); labs('$A$','$R,A~(\\alpha = 0.5)$\n\n$R$')
  sca(ah[1][1]); plotBb(MA,MR,MR,MA,Ma); Rticks(); Aticks(); labs('$A$',None)
  sca(ah[1][2]); plotRB(MA,MR,MR,MA,Ma); Rticks(); Aticks(); labs('$A$',None)
  # A vs a (R<1)
  MA,Ma = np.meshgrid(A,a); MR = np.log10(.1)
  sca(ah[2][0]); plotBw(MA,Ma,MR,MA,Ma); Aticks(); labs('$A$','$\\alpha,A~(R = 0.1)$\n\n$\\alpha$')
  sca(ah[2][1]); plotBb(MA,Ma,MR,MA,Ma); Aticks(); labs('$A$',None)
  sca(ah[2][2]); plotRB(MA,Ma,MR,MA,Ma); Aticks(); labs('$A$',None)
  # A vs a (R > 1)
  MA,Ma = np.meshgrid(A,a); MR = np.log10(5)
  sca(ah[3][0]); plotBw(MA,Ma,MR,MA,Ma); Aticks(); labs('$A$','$\\alpha,A~(R = 5)$\n\n$\\alpha$')
  sca(ah[3][1]); plotBb(MA,Ma,MR,MA,Ma); Aticks(); labs('$A$',None)
  sca(ah[3][2]); plotRB(MA,Ma,MR,MA,Ma); Aticks(); labs('$A$',None)
  # clean up
  plt.tight_layout()
  fh.subplots_adjust(bottom=.13)
  cbar('inferno',LB,orientation='horizontal',pos=[.12,.05,.54,.03],ticks=LB[::4],labels=10**LB[::4])
  cbar('viridis',LR,orientation='horizontal',pos=[.74,.05,.22,.03],ticks=LR[::2],labels=LR[::2])
  plt.savefig('B.xph.surf.pdf')
  # plt.show()

# maxBratio(*3*[2])
# maxBratio(*3*[.5])
maxplot()
surface()
