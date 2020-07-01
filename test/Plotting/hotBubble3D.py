
import numpy as np
import os
import matplotlib.pyplot as plt


plt.rcParams['font.size']=20

dirin='./hotBubble3D/'

dirout='./hotBubble3D/'

#time step
t=0
# plot choice (xy, xz, yz)
#choice = 'xz'


nx=128
ny=128
nz=128
nvar=7 #6 in Jan's setup the output was 6 in the new input files it is 7

# axis length
xmin = -4000
xmax =  4000
zmin =  0
zmax =  8000

start = 0
stop = 40
step = 1

ic = 0
for t in np.arange(start,stop,step):

  fL1 = open(dirin+'pf_all.dat', 'rb')
  fL1.seek(4*t*nx*ny*nz*nvar, os.SEEK_SET)
  f4d=np.fromfile(fL1, dtype='float32', count=nx*ny*nz*nvar)
  f4d=f4d.reshape((nvar,nz,ny,nx))
  fL1.close()
  
  dens = f4d[0,:,:,:]
  dens_pert = dens - np.mean(dens,2)
  u = f4d[1,:,:,:]
  v = f4d[2,:,:,:]
  w = f4d[3,:,:,:]
  exner = f4d[4,:,:,:]
  potTemp = f4d[5,:,:,:]
  potTemp_pert = potTemp - np.mean(potTemp,1)
  
  x = np.linspace(xmin, xmax, nx)
  #x=np.divide(x,1000)
  z=np.linspace(zmin, zmax, nz)
  #z=np.divide(z,1000)
  
  X, Z = np.meshgrid(x, z)
  
  print('u', np.amax(u[:,0,:]), np.amin(u[:,0,:]))
  print('v', np.amax(v[:,0,:]), np.amin(v[:,0,:]))
  print('w', np.amax(w[:,0,:]), np.amin(w[:,0,:]))
  print('dens_pert', np.amax(dens_pert[:,64,:]), np.amin(dens_pert[:,64,:]))
  print('exner', np.amax(exner[:,0,:]), np.amin(exner[:,0,:]))
  print('potTemp', np.amax(potTemp[:,0,:]), np.amin(potTemp[:,0,:]))
  print('potTemp_pert', np.amax(potTemp_pert[:,64,:]), np.amin(potTemp_pert[:,64,:]))
  
  dpi = 300
  #fig1 = plt.figure(figsize=(6600/dpi, 6600/dpi), dpi=dpi)
  fig1 = plt.figure()
  
  plt.set_cmap('bwr')
  
  plt.xlim([-4000,4000])
  plt.ylim([0,8000])
  plt.ylabel('z[m]')
  plt.xlabel('x[m]')
  clevs = np.arange(-0.02,0.006,0.001)
  cnplot = plt.contourf(X,Z,dens_pert[:,64,:],clevs)
  fig1.colorbar(cnplot)
  plt.title('Density Perturbation')
  #-- save figure
  if ic < 10:
     stradd = '00'
  elif ic < 100:
     stradd = '0'
  else:
     stradd = '' 
  plt.savefig(dirout+'dens_'+stradd+str(t)+'.png', bbox_inches='tight', dpi=dpi)
  ic = ic + 1















