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
data = tools.ModelOutput(data_path + 'coldBubble/')
reference = tools.ModelOutput(reference_path + 'coldBubble/')

# Set slices.
it = - 1
ix = slice(256, 448)

print(" ".join(("Time:", str(data.variables['t'][it]), "s")))

# Set grid.
xx = data.variables['x'][ix] * 0.001
zz = data.variables['z'][:] * 0.001

# Get potential temperature fluctuations.
thetap = data.groups['atmvar'].variables['thetap'][it, :, 0, ix]
thetapref = reference.groups['atmvar'].variables['thetap'][it, :, 0, ix]
deltatheta = thetap - thetapref

# Create plot.
(levels, colormap) = style.symmetric_contours(thetap.min(), thetap.max())
(figure, axes) = plt.subplots()
plot = axes.contourf(xx, zz, thetap, levels, cmap = colormap)
axes.set_xlabel(r"$x \, \mathrm{[km]}$")
axes.set_ylabel(r"$z \, \mathrm{[km]}$")
figure.colorbar(plot, label = r"$\theta' \, \mathrm{[K]}$")
figure.savefig("".join((data_path, "/results/coldBubble.png")), dpi = 500)

# Create difference plot.
if data_path != reference_path:
  (levels, colormap) = style.symmetric_contours(deltatheta.min(), \
      deltatheta.max())
  (figure, axes) = plt.subplots()
  plot = axes.contourf(xx, zz, deltatheta, levels, cmap = colormap)
  axes.set_xlabel(r"$x \, \mathrm{[km]}$")
  axes.set_ylabel(r"$z \, \mathrm{[km]}$")
  figure.colorbar(plot, label = r"$\Delta \theta' \," \
      r"\mathrm{\left[K\right]}$")
  figure.savefig("".join((data_path, "/results/coldBubble_difference.png")), \
      dpi = 500)