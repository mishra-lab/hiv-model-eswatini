import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib as mpl
np.set_printoptions(suppress=True)

# comparing how per-contact transmission modifiers are applied
# to compute per-partnership probability of transmission:
# - acsp: all contacts in some partnerships (average)
# - scap: some contacts in all partnerships (binomial)
# beta - per-contact probability
# A    - total contacts per partnership
# R    - relative beta for the factor
# a    - (\alpha) proportion of contacts affected by R
# maxBratio: finds conditions (R,A,a) to maximize Bscap / Bacsp
# plot: big grid of Bscap, Bacsp, and Bratio for pairs of R,A,a

beta = .0034 # Boily2009

def Bacsp(R,A,a):
  return (1-(1-R*beta)**A)*(a) + (1-(1-beta)**A)*(1-a)

def Bscap(R,A,a):
  return 1 - ((1-R*beta)**(A*(a)) * (1-beta)**(A*(1-a)))

def Bratio(R,A,a):
  return Bscap(R,A,a) / Bacsp(R,A,a)

def maxBratio(R0=.5):
  x0 = [R0,50,.5] # initial conditions
  xb = [(.01,10),(5,500),(.1,.9)] # bounds
  jfun = lambda x: -Bratio(*np.power(10,x)) # maximize ratio
  M = minimize(jfun,np.log10(x0),method='L-BFGS-B',bounds=np.log10(xb))
  xopt = np.power(10,M.x)
  print('f = {:.3f}\n  R = {:.2f}\n  A = {:0f}\n  a = {:3f}'.format(-M.fun,*xopt)+
        '\nB_acsp = {:.4f}\nB_scap = {:.4f}'.format(Bacsp(*xopt),Bscap(*xopt)))

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

def plot():
  Ni = 32; Nc = 10
  R = np.linspace( -2,  1,Ni)
  A = np.linspace(  3,300,Ni)
  a = np.linspace(  0,  1,Ni)
  fh,ah = plt.subplots(3,4,figsize=(13,9))
  Rticks = lambda: plt.yticks(ticks=[-2,-1,0,1],labels=[.01,.1,1,10])
  plotratio = lambda M1,M2,MR,MA,Ma: contour(M1,M2,Bratio(10**MR,MA,Ma),Nc,vmin=1,vmax=1.5,cmap='viridis')
  plotscap  = lambda M1,M2,MR,MA,Ma: contour(M1,M2,Bscap(10**MR,MA,Ma),Nc,vmin=0,vmax=1,cmap='inferno')
  plotacsp  = lambda M1,M2,MR,MA,Ma: contour(M1,M2,Bacsp(10**MR,MA,Ma),Nc,vmin=0,vmax=1,cmap='inferno')
  smr = plt.cm.ScalarMappable(cmap=plt.get_cmap('viridis',Nc),norm=mpl.colors.Normalize(vmin=1,vmax=1.5))
  smx = plt.cm.ScalarMappable(cmap=plt.get_cmap('inferno',Nc),norm=mpl.colors.Normalize(vmin=0,vmax=1))
  # R vs a
  Ma,MR = np.meshgrid(a,R); MA = np.median(A)
  sca(ah[0][0]); plotratio(Ma,MR,MR,MA,Ma); Rticks(); labs('$\\alpha$','Ratio: $B_1$ / $B_2$\n\n$R$','$R$ vs $\\alpha$')
  sca(ah[1][0]); plotscap (Ma,MR,MR,MA,Ma); Rticks(); labs('$\\alpha$','$B_1$: some contacts in all partnerships\n\n$R$')
  sca(ah[2][0]); plotacsp (Ma,MR,MR,MA,Ma); Rticks(); labs('$\\alpha$','$B_2$: all contacts in some partnerships\n\n$R$')
  # R vs A
  MA,MR = np.meshgrid(A,R); Ma = np.median(a)
  sca(ah[0][1]); plotratio(MA,MR,MR,MA,Ma); Rticks(); labs('$A$','$R$','$R$ vs $A$')
  sca(ah[1][1]); plotscap (MA,MR,MR,MA,Ma); Rticks(); labs('$A$','$R$')
  sca(ah[2][1]); plotacsp (MA,MR,MR,MA,Ma); Rticks(); labs('$A$','$R$')
  # A vs a (R<1)
  MA,Ma = np.meshgrid(A,a); MR = np.log10(.1)
  sca(ah[0][2]); plotratio(MA,Ma,MR,MA,Ma); labs('$A$','$\\alpha$','$\\alpha$ vs $A$ ($R < 1$)')
  sca(ah[1][2]); plotscap (MA,Ma,MR,MA,Ma); labs('$A$','$\\alpha$')
  sca(ah[2][2]); plotacsp (MA,Ma,MR,MA,Ma); labs('$A$','$\\alpha$')
  # A vs a (R > 1)
  MA,Ma = np.meshgrid(A,a); MR = np.log10(5)
  sca(ah[0][3]); plotratio(MA,Ma,MR,MA,Ma); labs('$A$','$\\alpha$','$\\alpha$ vs $A$ ($R > 1$)')
  sca(ah[1][3]); plotscap (MA,Ma,MR,MA,Ma); labs('$A$','$\\alpha$')
  sca(ah[2][3]); plotacsp (MA,Ma,MR,MA,Ma); labs('$A$','$\\alpha$')
  # clean up
  plt.tight_layout()
  fh.subplots_adjust(right=.87)
  plt.colorbar(smr,cax=fh.add_axes([.91,.72,.04,.23]))
  plt.colorbar(smx,cax=fh.add_axes([.91,.07,.04,.55]))
  plt.savefig('B.mod.vs.pdf')
  # plt.show()

maxBratio(R0=.5)
maxBratio(R0=2)
plot()