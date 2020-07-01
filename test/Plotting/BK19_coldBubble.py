import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import ndimage
import os

# LES-Plot

# grid
nx =1024#2048
ny =1
nz =128#256


# axis length
xmin = -25.6
xmax = 25.6
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


dens = f2d[5,0:200,nx/2:nx/2+800]
print(max(dens[0,:]), min(dens[0,:]))
print([ n for n,i in enumerate(dens[0,:]) if i<(-1.) ][-1])



x = np.linspace(0, 20., 800)
y = np.linspace(ymin, 5., 200)
X, Y = np.meshgrid(x, y)

print(x[[ n for n,i in enumerate(dens[0,:]) if i<(-1.) ][-1]])


#x = np.linspace(0,15000,300)

# plt.plot(x,dens[24,:])
# plt.show

ncont = 20
cmap = 'gist_earth'


plt.figure(figsize=(20,5))
plt.xlabel(xlab,size=15)
plt.ylabel(ylab,size=15)
#plt.xticks(np.arange(0.,xmax,step=2.))
#plt.title('Heating Term')
#cs4 = plt.contourf(X, Y, dens,ncont , cmap=cmap) 

cs4 = plt.contourf(X, Y, dens, np.arange(-16.5,1.5,1.), cmap=cmap)
plt.contour(cs4,2, colors='black')
cb = plt.colorbar(cs4)
cb.ax.tick_params(labelsize=15)
plt.xticks(size=15)
plt.yticks(size=15)
plt.show()

