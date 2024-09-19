import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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

x = data.variables['x'][:] * 0.001 - 30
zz = data.variables['z'][:, :, 0, :] * 0.001
xx, _ = np.meshgrid(x, zz[- 1, :, 0])

dz = zz[0, 1, 0] - zz[0, 0, 0]

nt = len(data.dimensions['time'])

# Print last output time.
print(" ".join(("Time:", str(data.variables['time'][- 1] / 60 / 60), "h")))

# Set vertical velocities.
ww = data.groups['atmvar'].variables['w'][:, :, 0, :]
wref = reference.groups['atmvar'].variables['w'][:, :, 0, :]

# Compute difference.
deltaw = ww - wref

# Compute maximum and contour levels.
peak = np.max(np.abs(ww[- 1, :, :]))
levels = np.arange(- peak, peak, 0.25)

# Make plot.
fig, axes = plt.subplots()
plot = axes.pcolormesh(xx, zz[- 1, :, :], ww[- 1, :, :], cmap = 'seismic', \
    vmin = - peak, vmax = peak, shading = 'gouraud')
axes.contour(xx, zz[- 1, :, :], ww[- 1, :, :], levels = levels, colors \
    = 'black', linewidths = 1.0)
axes.plot(xx[0, :], zz[- 1, 0, :] - dz / 2., color = 'black', linewidth = 1.0)
axes.set_xlim(- 10.0, 10.0)
axes.set_ylim(0.0, 10.0)
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
fig.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
fig.savefig("".join((data_path, "/results/mountainwave.pdf")))
fig.savefig("".join((data_path, "/results/mountainwave.png")), dpi = 500)

# Make animation.
if make_animation:

  fig, axes = plt.subplots()
  plot = axes.pcolormesh(xx, zz[- 1, :, :], ww[- 1, :, :], cmap = 'seismic', \
      vmin = - peak, vmax = peak, shading = 'gouraud')
  axes.contour(xx, zz[- 1, :, :], ww[- 1, :, :], levels = levels, colors \
      = 'black', linewidths = 1.0)
  axes.plot(xx[0, :], zz[- 1, 0, :] - dz / 2., color = 'black', linewidth \
      = 1.0)
  axes.set_xlim(- 10.0, 10.0)
  axes.set_ylim(0.0, 10.0)
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")

  def animate(it):
    axes.clear()
    axes.pcolormesh(xx, zz[it, :, :], ww[it, :, :], cmap = 'seismic', vmin = \
        - peak, vmax = peak, shading = 'gouraud')
    axes.contour(xx, zz[it, :, :], ww[it, :, :], levels = levels, colors \
        = 'black', linewidths = 1.0)
    axes.plot(xx[0, :], zz[it, 0, :] - dz / 2., color = 'black', linewidth \
        = 1.0)
    axes.set_xlim(- 10.0, 10.0)
    axes.set_ylim(0.0, 10.0)
    axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
    axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
    return axes

  fig.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
  animation_output = animation.FuncAnimation(fig, animate, frames = nt)
  writer = animation.FFMpegWriter(fps = 10)
  animation_output.save("".join((data_path, "/results/mountainwave.avi")), \
      writer = writer)

# Make difference plot.
if data_path != reference_path:
  peak = np.max(np.abs(deltaw[- 1, :, :]))
  fig, axes = plt.subplots()
  plot = axes.pcolormesh(xx, zz[- 1, :, :], deltaw[- 1, :, :], cmap \
      = 'seismic', vmin = - peak, vmax = peak, shading = 'gouraud')
  axes.contour(xx, zz[- 1, :, :], deltaw[- 1, :, :], colors = 'black', \
      linewidths = 1.0)
  axes.plot(xx[0, :], zz[- 1, 0, :] - dz / 2., color = 'black', linewidth \
      = 1.0)
  axes.set_xlim(- 10.0, 10.0)
  axes.set_ylim(0.0, 10.0)
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
  fig.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
  fig.savefig("".join((data_path, "/results/mountainwave_difference.pdf")))
  fig.savefig("".join((data_path, "/results/mountainwave_difference.png")), dpi = 500)