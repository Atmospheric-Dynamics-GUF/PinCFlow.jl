import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tools
#import style

from  data_path import *

print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + '/wavePacket3D')
reference = tools.ModelOutput(reference_path + '/wavePacket3D')

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
#(levels, colormap) = style.symmetric_contours(rhop.min(), rhop.max())
levels = np.linspace(rhop.min(), rhop.max())
colormap = "seismic"
(fig, axes) = plt.subplots()
plot = axes.contourf(xx, zz, rhop, levels, cmap = colormap)
axes.set_xlim(2000.0, 7000.0)
axes.set_ylim(10.0, 50.0)
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
fig.colorbar(plot, label = r"$\rho' \, \mathrm{\left[kg \, m^{- 3}\right]}$")
fig.savefig("../results/wavePacket3D.pdf")


# Create difference plot.
if (not (rhop == rhopref).all() ):
  (fig, axes) = plt.subplots()
  levels = np.linspace(deltarho.min(), deltarho.max())
  plot = axes.contourf(xx, zz, deltarho, levels, cmap = colormap)
  axes.set_xlim(2000.0, 7000.0)
  axes.set_ylim(10.0, 50.0)
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
  fig.colorbar(plot, label = r"$\rho' \, \mathrm{\left[kg \, m^{- 3}\right]}$")
  fig.savefig("../results/wavePacket3D_difference.pdf")
