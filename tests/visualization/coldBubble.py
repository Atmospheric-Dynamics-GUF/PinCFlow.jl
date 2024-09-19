import subprocess
import numpy as np
import matplotlib.pyplot as plt
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
        reference_path = "/scratch/atmodynamics/" + user_name + \
            "/PF/pinc/reference"

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

nx = len(data.dimensions['x'])
nz = len(data.dimensions['z'])

# Select area.
xx = data.variables['x'][int(0.5 * nx):int(0.875 * nx)] * 0.001 + 25
zz = data.variables['z'][:] * 0.001

theta = data.groups['atmvar'].variables['thetap'][it, :, 0, int(0.5 \
    * nx):int(0.875 * nx)]
thetaref = reference.groups['atmvar'].variables['thetap'][it, :, 0, int(0.5 \
    * nx):int(0.875 * nx)]

xlabel = r"$x \, \mathrm{[km]}$"
ylabel = r"$z \, \mathrm{[km]}$"

# Compute difference.
deltatheta = theta - thetaref

# Make plot.
maximum = np.max(np.abs(theta))
figure, axes = plt.subplots()
plot = axes.pcolormesh(xx, zz, theta, vmin = - maximum, vmax = maximum, \
    shading = "gouraud", cmap = "seismic")
axes.set_xlabel(xlabel)
axes.set_ylabel(ylabel)
figure.colorbar(plot, label = r"$\theta' \, \mathrm{[K]}$")
figure.savefig("".join((data_path, "/results/coldBubble.pdf")))
figure.savefig("".join((data_path, "/results/coldBubble.png")), dpi = 500)

# Make difference plot.
if data_path != reference_path:
  maximum = np.max(np.abs(deltatheta))
  figure, axes = plt.subplots()
  plot = axes.pcolormesh(xx, zz, deltatheta, vmin = - maximum, vmax = maximum, \
      shading = "gouraud", cmap = "seismic")
  axes.set_xlabel(xlabel)
  axes.set_ylabel(ylabel)
  figure.colorbar(plot, label = r"$\Delta \theta' \," r"\mathrm{\left[K\right" \
      r"]}$")
  figure.savefig("".join((data_path, "/results/coldBubble_difference.pdf")))
  figure.savefig("".join((data_path, "/results/coldBubble_difference.png")), \
      dpi = 500)
