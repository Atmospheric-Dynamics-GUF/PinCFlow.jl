import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tools
#import style

from  data_path import *
  
print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + '/mountainwave')
reference = tools.ModelOutput(reference_path + '/mountainwave')

# Set slices.
it = - 1

print(" ".join(("Time:", str(data.variables['t'][it]), "s")))

# Set grid.
xx = data.variables['x'][:] * 0.001 - 30
zz = data.variables['z'][it, :, 0] * 0.001
hh = data.variables['hm'][it, 0] * 0.001
xx = xx * np.ones_like(zz)

# Get vertical wind.
ww = data.groups['atmvar'].variables['w'][it, :, 0]
wref = reference.groups['atmvar'].variables['w'][it, :, 0]
deltaw = ww - wref

# Create plot.
#(levels, colormap) = style.symmetric_contours(ww.min(), ww.max())
levels = np.linspace(ww.min(), ww.max())
colormap = "seismic"
(fig, axes) = plt.subplots()
plot = axes.contourf(xx, zz, ww, levels, cmap = colormap)
axes.plot(xx[0], hh, color = 'black', linewidth = 1.0)
axes.set_xlim(- 10.0, 10.0)
axes.set_ylim(0.0, 10.0)
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
fig.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
fig.savefig("../results/mountainwave.pdf")

# Create difference plot.
if (not (ww == wref).all() ):
  maximum = numpy.max(numpy.abs(deltatheta))
  (fig, axes) = plt.subplots()
  plot = axes.contourf(xx, zz, deltaw, vmin = - maximum,
            vmax = maximum, shading = "gouraud", cmap = colormap)
  axes.plot(xx[0], hh, color = "black", linewidth = 1.0)
  axes.set_xlim(- 10.0, 10.0)
  axes.set_ylim(0.0, 10.0)
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
  fig.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
  fig.savefig("../results/mountainwave_difference.pdf")
