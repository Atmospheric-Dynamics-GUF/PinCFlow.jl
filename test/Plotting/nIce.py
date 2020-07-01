import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from PincFloitIO import PincFloitIO
import os
import sys
import matplotlib.colors as mcolors

###############################################################################

# parameters:
Mole_mass_water = 18.01528e-3 #kg/mol
Mole_mass_dryAir = 28.9644e-3 #kg/mol
epsilon0 = Mole_mass_water/Mole_mass_dryAir
p0 = 101325.0 #Pa
g = 9.81 #m/s^2
gamma = 5./3.
gamma_1 = gamma-1.
kappa = gamma_1/gamma
kappaInv = 1./kappa
Rsp = 287.0 #J/kg/K
T0 = 220 #K

# fitting
afit0=-62.192670609121556
afit1=254.77490427507394 
delta=1.522

###############################################################################
# functions:

def get_units(var):
  switcher = {
      'rho': r'$ \left[\frac{kg}{m^3}\right]$',
      'uu': r'$ \left[\frac{m}{s}\right]$',
      'vv': r'$ \left[\frac{m}{s}\right]$',
      'ww': r'$ \left[\frac{m}{s}\right]$',
      'pi': '',
      'theta': r'$ \left[K\right]$',
      'dynSma': '',
      'nAer': r'$ \left[\frac{1}{kg}\right]$',
      'nIce': r'$ \left[\frac{1}{kg}\right]$',
      'qIce': r'$ \left[\frac{kg}{kg}\right]$',
      'qv': r'$ \left[\frac{kg}{kg}\right]$',
      'SIce': '',
      'Temp': r'$ \left[K\right]$',
      'prs': r'$ \left[Pa\right]$',
      'mIce': r'$ \left[kg\right]$'
  }
  return switcher.get(var,"error")

def PStrat(k):
  return p0*np.exp(-data.zz[k]*g/(gamma*Rsp*T0)) #isothermal case

def pIce(T):
  return np.exp(9.550426-5723.265/T+3.53068*np.log(T)-0.00728332*T)

def p_saturation(T):
  return np.exp(54.842763-6763.22/T-4.210*np.log(T)+0.000367*T +np.tanh(0.0415*(T-218.8))*(53.878-1331.22/T-9.44523*np.log(T) +0.014025*T))

def SIce_threshold(T):
    awi = pIce(T) / p_saturation(T)
    return (-afit0 + delta)/(afit1 * awi) + 1.

def find_temperature():
  T = np.zeros((len_zz, len_xx))
  for i in range(len_xx):
    for k in range(len_zz):
      T[k,i] = (PStrat(k)/(Rsp*data.rho[timestep,k,0,i]) ) * ((PStrat(k)/p0)**gamma_1 + data.pi[timestep,k,0,i]) - T0
  return T

def find_SIce(): 
  SIce = np.zeros((len_zz, len_xx))
  for i in range(len_xx):
    for k in range(len_zz):
        T = (PStrat(k)/(Rsp*data.rho[timestep,k,0,i]) ) * ((PStrat(k)/p0)**gamma_1 +data.pi[timestep,k,0,i])
        p = p0 * ( (PStrat(k)/p0)**gamma_1  + data.pi[timestep,k,0,i] )**kappaInv
        SIce[k,i] = data.qv[timestep,k,0,i] * p / ( epsilon0 * pIce(T) ) #- SIce_threshold(T)
  return SIce

def find_prs():
  p = np.zeros((len_zz, len_xx)) 
  for i in range(len_xx): 
    for k in range(len_zz):
      p[k,i] = p0 * ((PStrat(k)/p0)**gamma_1 + data.pi[timestep,k,0,i])**kappaInv
  return p

###############################################################################

colors = [(0,0,0,c) for c in np.linspace(0,1,2)]
cmaptop = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=5)

pwd = os.popen("pwd").read()[0:-1]

data = PincFloitIO(pwd+'/pf_all.dat',pwd+'/input.f90',nvar=11,nframe=1)#3122)
var_list = ['rho','uu','vv','ww','pi','theta','dynSma','nAer','nIce','qIce','qv']
data.set_varnames(var_list)
len_tt = int(len(data.tt))
len_xx = int(len(data.xx))
len_zz = int(len(data.zz)*2./3.)
times = data.tt
topography = np.zeros((len_zz, len_xx))
for i in range(len_xx):
  for k in range(len_zz):
    if (data.uu[0,k,0,i]<0.01):
      topography[k,i] = 1.0

var_given = sys.argv[1]
timestep = int(sys.argv[2])
if (var_given in var_list): 
  var_list = [var_given]
elif (var_given == 'Temp'):
  var_list = ['Temp']
elif (var_given == 'SIce'):
  var_list = ['SIce']
elif (var_given == 'prs'):
  var_list = ['prs']
elif (var_given == 'log-nIce'):
  var_list = ['log-nIce']
elif (var_given == 'mIce'):
  var_list = ['mIce']
elif (var_given!="all"):
  print("Unknown variable name entered. Program terminated.")

for var_name in var_list:
 for i in range(1):
  timestep=timestep+i
  print("Creating "+var_name+"-plot")
  if (var_name=='Temp'):
    plot_var = find_temperature()
  elif (var_name=='SIce'):
    plot_var = find_SIce()
  elif (var_name=='prs'):
    plot_var = find_prs()
  elif (var_name=='log-nIce'):
    plot_var = np.log(1+data.nIce[timestep,0:len_zz,0,:])
  elif (var_name=='mIce'):
    plot_var = np.log(data.qIce[timestep,0:len_zz,0,:]/(data.nIce[timestep,0:len_zz,0,:]+1.e-20))
  else:
    plot_var = data.ret(var_name)[timestep,0:len_zz,0,:]
  var = plot_var[:,:]
  fig = plt.figure()
  scat = plt.pcolormesh(data.xx,data.zz[0:len_zz],var)
  plt.colorbar().set_label(var_name+get_units(var_name),rotation=90)
  quiv = plt.quiver(data.xx,data.zz[0:len_zz],data.uu[timestep,0:len_zz,0,:],data.ww[timestep,0:len_zz,0,:],color='black')
  plt.pcolormesh(data.xx,data.zz[0:len_zz],topography,cmap=cmaptop)
  plt.xlabel("x [m]")
  plt.ylabel("z [m]")
  plt.show()

