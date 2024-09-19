import subprocess
import numpy as np
import matplotlib.pyplot as plt
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

# Set plot parameters.
it = - 1

# Print time.
print(" ".join(("Time:", str(data.variables['time'][it]), "s")))

# Set fields of interest.
xx = data.variables['x'][:] * 0.001
zz = data.variables['z'][:] * 0.001
theta = data.groups['atmvar'].variables['thetap'][it, :, 0, :]
thetaref = reference.groups['atmvar'].variables['thetap'][it, :, 0, :]
deltatheta = theta - thetaref
xlabel = r"$x \, \mathrm{\left[km\right]}$"
ylabel = r"$z \, \mathrm{\left[km\right]}$"

# Make plot.
maximum = np.max(np.abs(theta))
fig, axes = plt.subplots()
plot = axes.pcolormesh(xx, zz, theta, vmin = - maximum, vmax = maximum, \
    shading = 'gouraud', cmap = 'seismic')
axes.set_xlabel(xlabel)
axes.set_ylabel(ylabel)
fig.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
fig.savefig("".join((data_path, "/results/IGW.pdf")))
fig.savefig("".join((data_path, "/results/IGW.png")), dpi = 500)

# Make difference plot.
if data_path != reference_path:
  maximum = np.max(np.abs(deltatheta))
  fig, axes = plt.subplots()
  plot = axes.pcolormesh(xx, zz, deltatheta, vmin = - maximum, vmax = maximum, \
      shading = 'gouraud', cmap = 'seismic')
  axes.set_xlabel(xlabel)
  axes.set_ylabel(ylabel)
  fig.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
  fig.savefig("".join((data_path, "/results/IGW_difference.pdf")))
  fig.savefig("".join((data_path, "/results/IGW_difference.png")), dpi = 500)
