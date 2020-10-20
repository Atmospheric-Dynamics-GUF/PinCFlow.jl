import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import ndimage
import os

# LES-Plot

res = 1
# grid
nx =int(1024/res)
ny =1
nz =int(128/res)


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
t= 0

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

f = open('pf_all_50m.dat','rb')

f.seek(4*t*nx*ny*nz*nvar, os.SEEK_SET)
f2d=np.fromfile(f, dtype='float32', count=nx*ny*nz*nvar)
f2d=f2d.reshape((nvar,nz,ny,nx))
f.close()

if choice == 'xy':
    f2d=f2d[:,kzplot,:,:]
    xlab='x [1000 km]'
    ylab='y [1000 km]'
elif choice == 'xz':
    f2d=f2d[:,:,jyplot,:]
    ny=nz
    ymin=zmin
    ymax=zmax
    xlab='x [km]'
    ylab='z [km]'
elif choice == 'yz':
    f2d = f2d[:,:,:,ixplot]
    nx=ny
    xmin=ymin
    xmax=ymax
    ny=nz
    ymin=zmin
    ymax=zmax
    xlab='y [1000 km]'
    ylab='z [km]'

dens = f2d[5,0:100/res,nx/2:nx/2+320/res]
print(np.amax(dens[:,:]), np.amin(dens[:,:]))
#print(dens[0,:])
#print([ n for n,i in enumerate(dens[0,:]) if i<=(-1.) ][-1])


x = np.linspace(0, 16, 320/res)
y = np.linspace(0, 5, 100/res)
X, Y = np.meshgrid(x, y)

#print(x[[ n for n,i in enumerate(dens[0,:]) if i<(-1.) ][-1]])


# x = np.linspace(0,15000,300)

# plt.plot(x,dens[24,:])
# plt.show

ncont = 20
cmap = 'YlGnBu_r'


plt.figure(figsize=(16,5))
plt.xlabel(xlab,size=30)
plt.ylabel(ylab,size=30)
#plt.title('semi-implicit, t = 0 sec', size=20)
#plt.xticks(np.arange(0.,xmax,step=2.))
#plt.title('Heating Term')
#cs4 = plt.contourf(X, Y, dens,ncont , cmap=cmap) 

cs4 = plt.contourf(X, Y, dens, np.arange(-16.5,1.,1.), cmap=cmap)
plt.contour(cs4,2, colors='black',linewidths=3)
cb = plt.colorbar(cs4)
cb.ax.tick_params(labelsize=30)
plt.xticks(size=30)
plt.yticks(size=30)
plt.ylim(0.,5.,1.)
plt.xlim(0.,16.,2.)
plt.tight_layout()
plt.show()

