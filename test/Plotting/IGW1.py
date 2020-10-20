import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from scipy import ndimage
import os


# LES-Plot

# grid
nx =300
ny =1
nz =10


# axis length
xmin = -150.#0.
xmax = 150.#300.
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

# for y-z plots:
ixplot = 32
# for x-z plots:
jyplot = 0
# for x-y plots:
kzplot = 14#48

#f2d=np.fromfile('pf_all.dat', dtype='float32').reshape((nt,nvar,nz,ny,nx))
#f2d=f2d[t,:,:,:,:]

f = open('pf_all.dat','rb')

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

dens = f2d[5,:,:]


x = np.linspace(xmin, xmax, nx)
y = np.linspace(ymin, ymax, ny)
X, Y = np.meshgrid(x, y)

ncont = 20
#cmap = cm.viridis



plt.figure(figsize=(15,5))
plt.xlabel(xlab,size=15)
plt.ylabel(ylab,size=15)
#plt.title('explicit, t = 3000sec',size = 20)
#plt.xticks(np.arange(0.,xmax,step=2.))
#plt.title('Heating Term')
#cs4 = plt.contourf(X, Y, dens,ncont , cmap=cmap) 

cs4 = plt.contourf(X, Y, dens,cmap = 'YlGnBu_r')
#cs4 = plt.contourf(X, Y, dens, np.arange(0.,0.0105,0.001), cmap ='YlGnBu_r')
cs4 = plt.contourf(X, Y, dens, np.arange(-0.003,0.0035,0.0005), cmap ='YlGnBu_r')
#YlGnBu_r
plt.contour(cs4,2, colors='black')
cb = plt.colorbar(cs4)
cb.ax.tick_params(labelsize=15)
plt.xticks(size=15)
plt.yticks(size=15)
plt.show()

