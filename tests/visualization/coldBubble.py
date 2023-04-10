import numpy
import matplotlib.pyplot as pyplot
import tools
import style

# Import data.
data = tools.ModelOutput("../coldBubble/")
reference = tools.ModelOutput("../coldBubble/")

# Set plot parameters.
choice = "xz"
it = - 1
ix = 16
iy = 0
iz = 7

# Print time.
print(" ".join(("Time:", str(data.tt[it]), "s")))

# Set fields of interest.
if choice == "xy":
    psi = data.psi[it, :, iz]
    psiref = reference.psi[it, :, iz]
    xx = 0.001 * data.xx[iz]
    yy = 0.001 * data.yy[iz]
    xlabel = r"$x \, \mathrm{\left[km\right]}$"
    ylabel = r"$y \, \mathrm{\left[km\right]}$"
elif choice == "xz":
    psi = data.psi[it, :, :, iy]
    psiref = reference.psi[it, :, :, iy]
    xx = 0.001 * data.xx[:, iy]
    yy = 0.001 * data.zz[:, iy]
    xlabel = r"$x \, \mathrm{\left[km\right]}$"
    ylabel = r"$z \, \mathrm{\left[km\right]}$"
elif choice == "yz":
    psi = data.psi[it, ..., ix]
    psiref = reference.psi[it, ..., ix]
    xx = 0.001 * data.yy[..., ix]
    yy = 0.001 * data.zz[..., ix]
    xlabel = r"$y \, \mathrm{\left[km\right]}$"
    ylabel = r"$z \, \mathrm{\left[km\right]}$"

# Select area.
xx = xx[:, int(0.5 * data.nx):int(0.875 * data.nx)]
yy = yy[:, int(0.5 * data.nx):int(0.875 * data.nx)]
psi = psi[..., int(0.5 * data.nx):int(0.875 * data.nx)]
psiref = psiref[..., int(0.5 * data.nx):int(0.875 * data.nx)]

# Compute difference.
deltapsi = psi - psiref

# Subtract offset.
psi[5] = psi[5] - 300.0

# Make plot.
maximum = numpy.max(numpy.abs(psi[5]))
figure, axes = pyplot.subplots()
plot = axes.pcolormesh(xx, yy, psi[5], vmin = - maximum, vmax = maximum,
        shading = "gouraud", cmap = "seismic")
axes.set_xlabel(xlabel)
axes.set_ylabel(ylabel)
figure.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
figure.savefig("../results/coldBubble.pdf")
figure.savefig("../results/coldBubble.png", dpi = 500)

# Make difference plot.
if (data.psi != reference.psi).all():
    maximum = numpy.max(numpy.abs(deltapsi[5]))
    figure, axes = pyplot.subplots()
    plot = axes.pcolormesh(xx, yy, deltapsi[5], vmin = - maximum,
            vmax = maximum, shading = "gouraud", cmap = "seismic")
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    figure.colorbar(plot, label = r"$\Delta \theta' \,"
            r"\mathrm{\left[K\right]}$")
    figure.savefig("../results/coldBubble_difference.pdf")
    figure.savefig("../results/coldBubble_difference.png", dpi = 500)
