import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tools
from scipy.ndimage import gaussian_filter
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
	reference_path = data_path

print("data_path =", data_path)
print("reference_path =", reference_path)

data = tools.ModelOutput(data_path + '/2DWP90_msgwam/')

xx = data.xx[:, 0, :] * 0.001
zz = data.zz[:, 0, :] * 0.001

# Large-scale tracer distribution change (tLS)
# and zonal, meridional, and vertical wind (uLS, vLS, wLS)
tLS = np.mean(data.psi[- 1, - 1, :, :, :], axis = 1)
uLS = np.mean(data.psi[- 1, 1, :, :, :], axis = 1)
vLS = np.mean(data.psi[- 1, 2, :, :, :], axis = 1)
wLS = np.mean(data.psi[- 1, 3, :, :, :], axis = 1)

print(" ".join(("Time:", str(np.round(data.tt[- 1] / 3600, 1)), "h")))

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
plt.colorbar(pplot00, ax = axes[0, 0], ticks = np.arange(- 0.044, 0.0441, \
    0.004)[1::2])

axes[0, 1].set_title(r'Large-scale meridional wind $v \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
pplot01 = axes[0, 1].contourf(xx, zz, vLS, cmap = 'seismic', vmin = - 0.02, \
    vmax = 0.02, levels = np.arange(- 0.02, 0.021, 0.002))
plt.colorbar(pplot01, ax = axes[0, 1], ticks = np.arange(- 0.02, 0.021, \
    0.002)[1::2])

axes[1, 0].set_title(r'Large-scale vertical wind $w \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
pplot10 = axes[1, 0].contourf(xx, zz, wLS, cmap = 'seismic', vmin = - 12e-5, \
    vmax = 12e-5, levels = np.arange(- 12e-5, 12.1e-5, 2e-5))
plt.colorbar(pplot10, ax = axes[1, 0], ticks = np.arange(- 12e-5, 12.1e-5, \
    2e-5)[1::2])

axes[1, 1].set_title('Large-scale tracer distribution')
pplot11 = axes[1, 1].contourf(xx, zz, tLS, cmap = 'seismic', vmin = - 1.2, \
    vmax = 1.2, levels = np.arange(- 1.2, 1.3, 0.1))
plt.colorbar(pplot11, ax = axes[1, 1], ticks = np.arange(- 1.2, 1.3, \
    0.1)[1::3])

plt.savefig("".join((data_path, "/results/2DWP90_msgwam.pdf")))
plt.savefig("".join((data_path, "/results/2DWP90_msgwam.png")), dpi = 500)

if data_path != reference_path:
	reference = tools.ModelOutput(reference_path + '/2DWP90_msgwam/')

	xxref = reference.xx[:, 0, :] * 0.001
	zzref = reference.zz[:, 0, :] * 0.001

	 # Large-scale tracer distribution change (tLS)
	 # and zonal, meridional, and vertical wind (uLS, vLS, wLS)
	tLSref = np.mean(reference.psi[- 1, - 1, :, :, :], axis = 1)
	uLSref = np.mean(reference.psi[- 1, 1, :, :, :], axis = 1)
	vLSref = np.mean(reference.psi[- 1, 2, :, :, :], axis = 1)
	wLSref = np.mean(reference.psi[- 1, 3, :, :, :], axis = 1)

	fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (8, 6))

	for ax in axes.flatten():
		ax.set_xlim(0, 9000)
		ax.set_ylim(15, 45)
		ax.set_xlabel('x [km]')
		ax.set_ylabel('z [km]')

	axes[0, 0].set_title(r'Large-scale zonal wind $u \, \left[\mathrm{m \,' \
      r's^{- 1}}\right]$')
	pplot00 = axes[0, 0].contourf(xx, zz, uLS - uLSref, cmap = 'seismic', \
      vmin = - 0.044, vmax = 0.044, levels = np.arange(- 0.044, 0.0441, 0.004))
	plt.colorbar(pplot00, ax = axes[0, 0], ticks = np.arange(- 0.044, 0.0441, \
    	0.004)[1::2])

	axes[0, 1].set_title(r'Large-scale meridional wind $v \, \left[\mathrm{m \,' \
      r's^{- 1}}\right]$')
	pplot01 = axes[0, 1].contourf(xx, zz, vLS - vLSref, cmap = 'seismic', \
      vmin = - 0.02, vmax = 0.02, levels = np.arange(- 0.02, 0.021, 0.002))
	plt.colorbar(pplot01, ax = axes[0, 1], ticks = np.arange(- 0.02, 0.021, \
    	0.002)[1::2])

	axes[1, 0].set_title(r'Large-scale vertical wind $w \, \left[\mathrm{m \,' \
      r's^{- 1}}\right]$')
	pplot10 = axes[1, 0].contourf(xx, zz, wLS - wLSref, cmap = 'seismic', \
      vmin = - 12e-5, vmax = 12e-5, levels = np.arange(- 12e-5, 12.1e-5, 2e-5))
	plt.colorbar(pplot10, ax = axes[1, 0], ticks = np.arange(- 12e-5, \
    	12.1e-5, 2e-5)[1::2])

	axes[1, 1].set_title('Large-scale tracer distribution')
	pplot11 = axes[1, 1].contourf(xx, zz, tLS - tLSref, cmap = 'seismic', \
      vmin = - 1.2, vmax = 1.2, levels = np.arange(- 1.2, 1.3, 0.1))
	plt.colorbar(pplot11, ax = axes[1, 1], ticks = np.arange(- 1.2, 1.3, \
    	0.1)[1::3])

	plt.savefig("".join((data_path, "/results/2DWP90_msgwam.pdf")))
	plt.savefig("".join((data_path, "/results/2DWP90_msgwam.png")), dpi = 500)
