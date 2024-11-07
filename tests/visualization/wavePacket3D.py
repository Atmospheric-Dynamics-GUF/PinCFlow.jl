import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tools
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
  data_path = "../"
  reference_path = "../"

print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + 'wavePacket3D/')
reference = tools.ModelOutput(reference_path + 'wavePacket3D/')

# Set slices.
it = - 1

print(" ".join(("Time:", str(data.variables['t'][it]), "s")))

# Set grid.
xx = data.variables['x'][:] * 0.001
zz = data.variables['z'][:] * 0.001

# Get density fluctuations.
rhop = data.groups['atmvar'].variables['rhop'][it, :, 0]
rhopref = reference.groups['atmvar'].variables['rhop'][it, :, 0]
deltarho = rhop - rhopref

# Create plot.
(levels, colormap) = style.symmetric_contours(rhop.min(), rhop.max())
(fig, axes) = plt.subplots()
plot = axes.contourf(xx, zz, rhop, levels, cmap = colormap)
axes.set_xlim(2000.0, 7000.0)
axes.set_ylim(10.0, 50.0)
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
fig.colorbar(plot, label = r"$\rho' \, \mathrm{\left[kg \, m^{- 3}\right]}$")
fig.savefig("".join((data_path, "/results/wavePacket3D.pdf")))
fig.savefig("".join((data_path, "/results/wavePacket3D.png")), dpi = 500)

# Create difference plot.
if data_path != reference_path:
  (levels, colormap) = style.symmetric_contours(deltarho.min(), deltarho.max())
  (fig, axes) = plt.subplots()
  plot = axes.contourf(xx, zz, deltarho, levels, cmap = colormap)
  axes.set_xlim(2000.0, 7000.0)
  axes.set_ylim(10.0, 50.0)
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
  fig.colorbar(plot, label = r"$\rho' \, \mathrm{\left[kg \, m^{- 3}\right]}$")
  fig.savefig("".join((data_path, "/results/wavePacket3D_difference.pdf")))
  fig.savefig("".join((data_path, "/results/wavePacket3D_difference.png")), \
      dpi = 500)