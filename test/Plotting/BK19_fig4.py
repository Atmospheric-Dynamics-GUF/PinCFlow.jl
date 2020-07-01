import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import ndimage
import os
import sys

# LES-Plot

# grid
nx =2048
ny =1
nz =256


# axis length
xmin = -25.6
xmax = 25.6
ymin = 0.0
ymax = 50.0
zmin = 0.
zmax = 6.4

# plot range
x_range = (0.0, 1.0)
y_range = (0.0, 1.0)
z_range = (0.0, 1.0)

# time step (beginning from 0!)
t= 1

# number of variables
nvar = 7

# plot choice (xy, xz, yz)
choice = 'xz'

# for y-z plots:
ixplot = 32
# for x-z plots:
jyplot = 0
# for x-y plots:
kzplot = 14#48

#f2d=np.fromfile('pf_all.dat', dtype='float32').reshape((nt,nvar,nz,ny,nx))
#f2d=f2d[t,:,:,:,:]

f400 = open('pf_all_400m.dat','rb')
f200 = open('pf_all_200m.dat','rb')
f100 = open('pf_all_100m.dat','rb')
f50 = open('pf_all_50m.dat','rb')
f25 = open('pf_all_25m.dat','rb')

f400.seek(4*t*nx/16*ny*nz/16*nvar, os.SEEK_SET)
f2d400=np.fromfile(f400, dtype='float32', count=nx/16*ny*nz/16*nvar)
f2d400=f2d400.reshape((nvar,nz/16,ny,nx/16))
f400.close()

f200.seek(4*t*nx/8*ny*nz/8*nvar, os.SEEK_SET)
f2d200=np.fromfile(f200, dtype='float32', count=nx/8*ny*nz/8*nvar)
f2d200=f2d200.reshape((nvar,nz/8,ny,nx/8))
f200.close()

f100.seek(4*t*nx/4*ny*nz/4*nvar, os.SEEK_SET)
f2d100=np.fromfile(f100, dtype='float32', count=nx/4*ny*nz/4*nvar)
f2d100=f2d100.reshape((nvar,nz/4,ny,nx/4))
f100.close()

f50.seek(4*t*nx/2*ny*nz/2*nvar, os.SEEK_SET)
f2d50=np.fromfile(f50, dtype='float32', count=nx/2*ny*nz/2*nvar)
f2d50=f2d50.reshape((nvar,nz/2,ny,nx/2))
f50.close()

f25.seek(4*t*nx*ny*nz*nvar, os.SEEK_SET)
f2d25=np.fromfile(f25, dtype='float32', count=nx*ny*nz*nvar)
f2d25=f2d25.reshape((nvar,nz,ny,nx))
f25.close()

f2d400=f2d400[:,:,jyplot,:]
f2d200=f2d200[:,:,jyplot,:]
f2d100=f2d100[:,:,jyplot,:]
f2d50=f2d50[:,:,jyplot,:]
f2d25=f2d25[:,:,jyplot,:]
ny=nz
ymin=zmin
ymax=zmax
xlab='x [km]'
ylab='z [km]'


dens400 = f2d400[5,:,nx/(2*16):nx/(2*16)+37]
dens200 = f2d200[5,:,nx/(2*8):nx/(2*8)+75]
dens100 = f2d100[5,:,nx/(2*4):nx/(2*4)+150]
dens50 = f2d50[5,:,nx/(2*2):nx/(2*2)+300]
dens25 = f2d25[5,:,nx/2:nx/2+600]

x400 = np.linspace(0, 15.,37.5 )
x200 = np.linspace(0, 15.,75 )
x100 = np.linspace(0, 15., 150)
x50 = np.linspace(0, 15., 300)
x25 = np.linspace(0, 15., 600)


plt.plot(x400,dens400[3,:], label = '400 m',color='black',linewidth=1)
plt.plot(x200,dens200[6,:], label = '200 m',color='red',ls='dashed',linewidth=3)
plt.plot(x100,dens100[12,:], label = '100 m',color='blue', ls ='dashdot',linewidth=3)
plt.plot(x50,dens50[24,:], label = '50 m',color='m',ls='solid',marker='o',linewidth=3)
plt.plot(x25,dens25[48,:], label = '25 m',color='green',ls='solid',marker='x',linewidth=2)
plt.legend(loc='lower left')
plt.xlabel('x[km]',size = 15)
plt.ylabel('theta - theta0 [K]',size=15)
#plt.ylim(-6.,1.,1.)
plt.xlim(0.,16.,5.)
plt.xticks(size=15)
plt.yticks(size=15)
plt.show()

