import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tools
#import style

from  data_path import *

print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + '/coldBubble')
reference = tools.ModelOutput(reference_path + '/coldBubble')

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
#(levels, colormap) = style.symmetric_contours(thetap.min(), thetap.max())
levels = np.linspace(thetap.min(), thetap.max())
colormap = "seismic"
(figure, axes) = plt.subplots()
plot = axes.contourf(xx, zz, thetap, levels, cmap = colormap)
axes.set_xlabel(r"$x \, \mathrm{[km]}$")
axes.set_ylabel(r"$z \, \mathrm{[km]}$")
figure.colorbar(plot, label = r"$\theta' \, \mathrm{[K]}$")
figure.savefig("../results/coldBubble.pdf")

# Create difference plot.
if (not (thetap == thetapref).all() ):
  print('there is some difference !!!')
  (figure, axes) = plt.subplots()
  levels = np.linspace(deltatheta.min(), deltatheta.max())
  plot = axes.contourf(xx, zz, deltatheta, levels, cmap = colormap)
  axes.set_xlabel(r"$x \, \mathrm{[km]}$")
  axes.set_ylabel(r"$z \, \mathrm{[km]}$")
  figure.colorbar(plot, label = r"$\Delta \theta' \," \
      r"\mathrm{\left[K\right]}$")
  figure.savefig("../results/coldBubble_difference.pdf")

