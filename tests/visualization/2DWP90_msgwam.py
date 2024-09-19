import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import netCDF4 as nc
import style

# Get host and user name.
host_name = subprocess.getoutput("hostname")
user_name = subprocess.getoutput("whoami")

if "levante" in host_name:
	 # Levante cluster
	data_path = "/scratch/b/" + user_name + "/PF/runs"
	reference_path = "/scratch/b/" + user_name + "/PF/pinc/reference"

elif "login" in host_name:
	 # Goethe cluster
	data_path = "/scratch/atmodynamics/" + user_name + "/PF/runs"
	reference_path = "/scratch/atmodynamics/" + user_name + "/PF/pinc/reference"

else:
	 # Local machine
	data_path = ".."
	reference_path = ".."

print("data_path =", data_path)
print("reference_path =", reference_path)

data = nc.Dataset(data_path + '/pincflow_data.nc')
reference = nc.Dataset(reference_path + '/pincflow_data.nc')

# large-scale tracer distribution change (tLS) 
# and zonal, meridional, and vertical wind (uLS, vLS, wLS)
tLS = data.groups['atmvar'].variables['tmrd'][-1, :, 0, :]
uLS = data.groups['atmvar'].variables['u'][-1, :, 0, :]
vLS = data.groups['atmvar'].variables['v'][-1, :, 0, :]
wLS = data.groups['atmvar'].variables['w'][-1, :, 0, :]

tLSref = reference.groups['atmvar'].variables['tmrd'][-1, :, 0, :]
uLSref = reference.groups['atmvar'].variables['u'][-1, :, 0, :]
vLSref = reference.groups['atmvar'].variables['v'][-1, :, 0, :]
wLSref = reference.groups['atmvar'].variables['w'][-1, :, 0, :]

deltat = tLS - tLSref 
deltau = uLS - uLSref 
deltav = vLS - vLSref 
deltaw = wLS - wLSref

x = data.variables['x'][:] * 0.001
z = data.variables['z'][:] * 0.001

xx, zz = np.meshgrid(x, z)


print(" ".join(("Time:", str(np.round(data.variables['time'][-1] / 3600, 1)), "h")))

fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (8, 6))

for ax in axes.flatten():
	ax.set_xlim(0, 9000)
	ax.set_ylim(15, 45)
	ax.set_xlabel('x [km]')
	ax.set_ylabel('z [km]')

axes[0, 0].set_title(r'Large-scale zonal wind $u \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
pplot00 = axes[0, 0].contourf(xx, zz, uLS, cmap = 'seismic', vmin = - 0.044, \
    vmax = 0.044, levels = np.arange(- 0.044, 0.0441, 0.004))
fig.colorbar(pplot00, ax = axes[0, 0], ticks = np.arange(- 0.044, 0.0441, \
    0.004)[1::2])

axes[0, 1].set_title(r'Large-scale meridional wind $v \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
pplot01 = axes[0, 1].contourf(xx, zz, vLS, cmap = 'seismic', vmin = - 0.02, \
    vmax = 0.02, levels = np.arange(- 0.02, 0.021, 0.002))
fig.colorbar(pplot01, ax = axes[0, 1], ticks = np.arange(- 0.02, 0.021, \
    0.002)[1::2])

axes[1, 0].set_title(r'Large-scale vertical wind $w \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
pplot10 = axes[1, 0].contourf(xx, zz, wLS, cmap = 'seismic', vmin = - 12e-5, \
    vmax = 12e-5, levels = np.arange(- 12e-5, 12.1e-5, 2e-5))
fig.colorbar(pplot10, ax = axes[1, 0], ticks = np.arange(- 12e-5, 12.1e-5, \
    2e-5)[1::2])

axes[1, 1].set_title('Large-scale tracer distribution')
pplot11 = axes[1, 1].contourf(xx, zz, tLS, cmap = 'seismic', vmin = - 1.2, \
    vmax = 1.2, levels = np.arange(- 1.2, 1.3, 0.1))
fig.colorbar(pplot11, ax = axes[1, 1], ticks = np.arange(- 1.2, 1.3, \
    0.1)[1::3])

fig.tight_layout()
fig.savefig("".join((data_path, "/results/2DWP90_msgwam.pdf")))
fig.savefig("".join((data_path, "/results/2DWP90_msgwam.png")), dpi = 500)


# Make difference plot.
if data_path != reference_path:
	peaku = np.max(np.abs(deltau))
	peakv = np.max(np.abs(deltav))
	peakw = np.max(np.abs(deltaw))
	peakt = np.max(np.abs(deltat))
	fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (8, 6))

	for ax in axes.flatten():
		ax.set_xlim(0, 9000)
		ax.set_ylim(15, 45)
		ax.set_xlabel('x [km]')
		ax.set_ylabel('z [km]')

	axes[0, 0].set_title(r'Large-scale zonal wind diff. $\Delta u \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
	pplot00 = axes[0, 0].contourf(xx, zz, deltau, cmap = 'seismic', 
							   vmin = - peaku, vmax = peaku)
	fig.colorbar(pplot00, ax = axes[0, 0])

	axes[0, 1].set_title(r'Large-scale meridional wind diff. $\Delta v \, \left[\mathrm{m \,' \
		r's^{- 1}}\right]$')
	pplot01 = axes[0, 1].contourf(xx, zz, deltav, cmap = 'seismic', 
							   vmin = - peakv, vmax = peakv)
	fig.colorbar(pplot01, ax = axes[0, 1])

	axes[1, 0].set_title(r'Large-scale vertical wind diff. $\Delta w \, \left[\mathrm{m \,' \
		r's^{- 1}}\right]$')
	pplot10 = axes[1, 0].contourf(xx, zz, deltaw, cmap = 'seismic', 
							   vmin = - peakw, vmax = peakw)
	fig.colorbar(pplot10, ax = axes[1, 0])

	axes[1, 1].set_title('Large-scale tracer distribution diff.')
	pplot11 = axes[1, 1].contourf(xx, zz, deltat, cmap = 'seismic', 
							   vmin = - peakt, vmax = peakt)
	fig.colorbar(pplot11, ax = axes[1, 1])

	fig.tight_layout()
	fig.savefig("".join((data_path, "/results/2DWP90_msgwam_difference.pdf")))
	fig.savefig("".join((data_path, "/results/2DWP90_msgwam_difference.png")), dpi = 500)
