
import numpy as np
import os
import matplotlib.pyplot as plt


plt.rcParams['font.size']=20

dirin='./1/'

dirout='./1/'

#time step
t=1
# plot choice (xy, xz, yz)
#choice = 'xz'


nx=512
ny=16
nz=1000
nvar=7 #6 in Jan's setup the output was 6 in the new input files it is 7

nx_wkb=128      
ny_wkb=1
nz_wkb=128
nvar_wkb=7 

# axis length
xmin = 0
xmax = 9000000
zmin = 0
zmax = 100000


fL1 = open(dirin+'pf_all.dat', 'rb')
fL1.seek(4*t*nx*ny*nz*nvar, os.SEEK_SET)
f4d=np.fromfile(fL1, dtype='float32', count=nx*ny*nz*nvar)
f4d=f4d.reshape((nvar,nz,ny,nx))
fL1.close()

dens = f4d[0,:,:,:]
u = f4d[1,:,:,:]
v = f4d[2,:,:,:]
w = f4d[3,:,:,:]
exner = f4d[4,:,:,:]
potTemp = f4d[5,:,:,:]

x = np.linspace(xmin, xmax, nx)
x=np.divide(x,1000)
z=np.linspace(zmin, zmax, nz)
z=np.divide(z,1000)

X, Z = np.meshgrid(x, z)

print('u', np.amax(u[:,0,:]), np.amin(u[:,0,:]))
print('v', np.amax(v[:,0,:]), np.amin(v[:,0,:]))
print('w', np.amax(w[:,0,:]), np.amin(w[:,0,:]))
print('dens', np.amax(dens[:,0,:]), np.amin(dens[:,0,:]))
print('exner', np.amax(exner[:,0,:]), np.amin(exner[:,0,:]))
print('potTemp', np.amax(potTemp[:,0,:]), np.amin(potTemp[:,0,:]))

dpi = 300
fig1 = plt.figure(figsize=(6600/dpi, 3300/dpi), dpi=dpi)

plt.set_cmap('bwr')

ax1 = fig1.add_subplot(231)
ax1.set_xlim([0,9000])
ax1.set_ylim([0,100])
ax1.set_ylabel('z[km]')
clevs1 = np.arange(-1,1,0.1)
cnplot1 = ax1.contourf(X,Z,u[:,0,:],clevs1)
fig1.colorbar(cnplot1)
ax1.set_title('U wind')

ax2 = plt.subplot(232)
ax2.set_xlim([0,9000])
ax2.set_ylim([0,100])
ax2.set_ylabel('z[km]')
clevs2 = np.arange(-1,1,0.1)
cnplot2 = ax2.contourf(X,Z,v[:,0,:],clevs2)
fig1.colorbar(cnplot2)
ax2.set_title('V wind')

ax3 = plt.subplot(233)
ax3.set_xlim([0,9000])
ax3.set_ylim([0,100])
ax3.set_ylabel('z[km]')
clevs3 = np.arange(-0.01,0.01,0.001)
cnplot3 = ax3.contourf(X,Z,w[:,0,:],clevs3)
fig1.colorbar(cnplot3)
ax3.set_title('W wind')

ax4 = plt.subplot(234)
ax4.set_xlim([0,9000])
ax4.set_ylim([0,100])
ax4.set_ylabel('z[km]')
ax4.set_xlabel('x[km]')
clevs4 = np.arange(-0.0002,0.0002,0.00005)
cnplot4 = ax4.contourf(X,Z,dens[:,0,:],clevs4)
fig1.colorbar(cnplot4)
ax4.set_title('Density')

ax5 = plt.subplot(235)
ax5.set_xlim([0,9000])
ax5.set_ylim([0,100])
ax5.set_ylabel('z[km]')
ax5.set_xlabel('x[km]')
clevs5 = np.arange(-0.01,0.01,0.001)
cnplot5 = ax5.contourf(X,Z,exner[:,0,:],clevs5)
fig1.colorbar(cnplot5)
ax5.set_title('Exner Pressure')

ax6 = plt.subplot(236)
ax6.set_xlim([0,9000])
ax6.set_ylim([0,100])
ax6.set_ylabel('z[km]')
ax6.set_xlabel('x[km]')
clevs6 = np.arange(300,3000,10)
cnplot6 = ax6.contourf(X,Z,potTemp[:,0,:],clevs6)
cbar = fig1.colorbar(cnplot6)
ax6.set_title('Potential Temperature')

plt.savefig(dirout+'wavePacket3D.png', bbox_inches='tight', dpi=dpi)















