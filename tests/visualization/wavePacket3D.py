import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import netCDF4 as nc
import style

# Set script parameter.
make_animation = False

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

# Import data.
data = nc.Dataset(data_path + '/pincflow_data.nc')
reference = nc.Dataset(reference_path + '/pincflow_data.nc')


it = -1 

# Print time.
print(" ".join(("Time:", str(data.variables['time'][it]), "s")))

# Set variables
x = data.variables['x'][:] * 0.001
z = data.variables['z'][:] * 0.001
xx, zz = np.meshgrid(x, z)
rhop = data.groups['atmvar'].variables['rhop'][-1, :, 0, :]
rhopref = reference.groups['atmvar'].variables['rhop'][-1, :, 0, :]

deltarho = rhop - rhopref 

rhopmax = np.max(np.abs(rhop))
deltamax = np.max(np.abs(deltarho))

fig, axes = plt.subplots()
plot = axes.pcolormesh(xx, zz, rhop, vmin = -rhopmax, vmax = rhopmax, 
					   shading = 'gouraud', cmap = 'seismic')
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
fig.colorbar(plot, label = r"$\rho' \," \
      r"\mathrm{\left[kg \, m^{- 3}\right]}$")
fig.savefig("".join((data_path, "/results/wavePacket3D.pdf")))
fig.savefig("".join((data_path, "/results/wavePacket3D.png")), dpi = 500)

if data_path != reference_path:
	fig, axes = plt.subplots()
	plot = axes.pcolormesh(xx, zz, deltarho, vmin = -deltamax, vmax = deltamax, 
						shading = 'gouraud', cmap = 'seismic')
	axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
	axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
	fig.colorbar(plot, label = r"$\rho' \," \
		r"\mathrm{\left[kg \, m^{- 3}\right]}$")
	fig.savefig("".join((data_path, "/results/wavePacket3D_difference.pdf")))
	fig.savefig("".join((data_path, "/results/wavePacket3D_difference.png")), dpi = 500)