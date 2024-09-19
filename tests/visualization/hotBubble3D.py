import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import netCDF4 as nc
import style

# Set script parameter.
make_animation = True

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

xx = data.variables['x'][:] * 0.001 + 4
zz = data.variables['z'][:] * 0.001

nt = len(data.dimensions['time'])
it = int(0.4 * len(data.dimensions['time']))
iy = int(0.5 * len(data.dimensions['y']))

# Print time.
print(" ".join(("Time:", str(data.variables['time'][it]), "s")))

# Set potential temperature.
theta = data.groups['atmvar'].variables['thetap'][:, :, iy, :]
thetaref = reference.groups['atmvar'].variables['thetap'][:, :, iy, :]

# Compute difference.
deltatheta = theta - thetaref

# Make plot.
fig, axes = plt.subplots()
peak = np.max(np.abs(theta[it, :, :]))
plot = axes.pcolormesh(xx, zz, theta[it, :, :], cmap = 'seismic', vmin = \
    - peak, vmax = peak, shading = 'gouraud')
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")

fig.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
fig.savefig("".join((data_path, "/results/hotBubble3D.pdf")))
fig.savefig("".join((data_path, "/results/hotBubble3D.png")), dpi = 500)

# Make difference plot.
if data_path != reference_path:
  fig, axes = plt.subplots()
  peak = np.max(np.abs(deltatheta[it, :, :]))
  plot = axes.pcolormesh(xx, zz, deltatheta[it, :, :], cmap = 'seismic', vmin \
      = - peak, vmax = peak, shading = 'gouraud')
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")

  fig.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
  fig.savefig("".join((data_path, "/results/hotBubble3D_difference.pdf")))
  fig.savefig("".join((data_path, "/results/hotBubble3D_difference.png")), dpi \
      = 500)

if make_animation:
  fig, axes = plt.subplots()
  peak = np.max(np.abs(theta[it, :, :]))
  plot = axes.pcolormesh(xx, zz, theta[0, :, :], cmap = 'seismic', vmin = \
      - peak, vmax = peak, shading = 'gouraud')
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")

  fig.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")

  def animate(nn):
    axes.clear()
    peak = np.max(np.abs(theta[it, :, :]))
    plot = axes.pcolormesh(xx, zz, theta[nn, :, :], cmap = 'seismic', vmin = \
        - peak, vmax = peak, shading = 'gouraud')
    axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
    axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
    return axes

  animation_output = animation.FuncAnimation(fig, animate, frames = nt)
  writer = animation.FFMpegWriter(fps = 10)
  animation_output.save("".join((data_path, "/results/hotBubble3D.avi")), \
      writer = writer)