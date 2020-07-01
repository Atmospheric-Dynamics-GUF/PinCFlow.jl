import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import ndimage
import os

# LES-Plot

# grid
nx =300
ny =1
nz =10


# axis length
xmin = 0.
xmax = 480000.
ymin = 0.0
ymax = 50.0
zmin = 0.
zmax = 10.

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

# for x-z plots:
jyplot = 0

#f2d=np.fromfile('pf_all.dat', dtype='float32').reshape((nt,nvar,nz,ny,nx))
#f2d=f2d[t,:,:,:,:]

f = open('pf_all.dat','rb')

f.seek(4*t*nx*ny*nz*nvar, os.SEEK_SET)
f2d=np.fromfile(f, dtype='float32', count=nx*ny*nz*nvar)
f2d=f2d.reshape((nvar,nz,ny,nx))
f.close()


f2d=f2d[:,:,jyplot,:]
ny=nz
ymin=zmin
ymax=zmax
xlab='x [km]'
ylab='z [km]'


dens = f2d[5,:,:]


x = np.linspace(xmin, xmax, nx)
y = np.linspace(ymin, ymax, ny)
X, Y = np.meshgrid(x, y)

ncont = 20
cmap = 'gist_earth'


plt.figure(figsize=(15,5))
plt.xlabel(xlab,size=15)
plt.ylabel(ylab,size=15)
#plt.xticks(np.arange(0.,xmax,step=2.))
#plt.title('Heating Term')
cs4 = plt.contourf(X, Y, dens,ncont , cmap=cmap) 

#cs4 = plt.contourf(X, Y, dens, np.arange(0.,0.011,0.001), cmap=cmap)
#cs4 = plt.contourf(X, Y, dens, np.arange(-0.0015,0.0035,0.0005), cmap=cmap)
plt.contour(cs4,2, colors='black')
cb = plt.colorbar(cs4)
cb.ax.tick_params(labelsize=15)
plt.xticks(size=15)
plt.yticks(size=15)
plt.show()

