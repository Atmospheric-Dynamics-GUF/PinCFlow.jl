import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tools
#import style

from  data_path import *
  
print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + '/IGW')
reference = tools.ModelOutput(reference_path + '/IGW')

# Set slices.
it = - 1

print(" ".join(("Time:", str(data.variables['t'][it]), "s")))

# Set grid.
xx = data.variables['x'][:] * 0.001
zz = data.variables['z'][:] * 0.001

# Get potential temperature fluctuations.
thetap = data.groups['atmvar'].variables['thetap'][it, :, 0]
thetapref = reference.groups['atmvar'].variables['thetap'][it, :, 0]
deltatheta = thetap - thetapref

# Create plot.
#(levels, colormap) = style.symmetric_contours(thetap.min(), thetap.max())
levels = np.linspace(thetap.min(), thetap.max())
colormap = "seismic"
(fig, axes) = plt.subplots()
plot = axes.contourf(xx, zz, thetap, levels, cmap = colormap)
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
fig.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
fig.savefig("../results/IGW.pdf")

# Create difference plot.
if (not (thetap == thetapref).all() ):
    print('there is some difference !!!')
    maximum = numpy.max(numpy.abs(deltatheta))
    figure, axes = pyplot.subplots()
    plot = axes.pcolormesh(xx, zz, deltatheta, vmin = - maximum,
            vmax = maximum, shading = "gouraud", cmap = "seismic")
    axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
    axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
    figure.colorbar(plot, label = r"$\Delta \theta' \,"
            r"\mathrm{\left[K\right]}$")
    fig.savefig("../results/IGW_difference.pdf")

