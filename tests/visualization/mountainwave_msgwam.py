import numpy
import matplotlib.pyplot as pyplot
import tools
#import style

data_path = "/scratch/atmodynamics/dolaptchiev/PF/runs"
ref_path = "/scratch/atmodynamics/dolaptchiev/PF/pinc/reference"

# Import data.
data = tools.ModelOutput(data_path+"/mountainwave_msgwam/")
reference = tools.ModelOutput(ref_path+"/mountainwave_msgwam/")

# Adust coordinate unit.
data.xx = 0.001 * data.xx
data.zz = 0.001 * data.zz

# Set plot parameters.
plot_all = False
it = 1

# Print time.
print(" ".join(("Time:", str(data.tt[it]), "s")))

# Compute difference.
deltapsi = data.psi - reference.psi

# Make plot.
if plot_all:
    umax = numpy.max(numpy.absolute(data.psi[it, 1, :, 0]))
    vmax = numpy.max(numpy.absolute(data.psi[it, 2, :, 0]))
    wmax = numpy.max(numpy.absolute(data.psi[it, 3, :, 0]))
    rhopmax = numpy.max(numpy.absolute(data.psi[it, 0, :, 0]))
    pipmax = numpy.max(numpy.absolute(data.psi[it, 4, :, 0]))
    figure, axes = pyplot.subplots(2, 3, figsize = (8.0, 6.0))
    zonal_wind = axes[0, 0].pcolormesh(data.xx[:, 0], data.zz[:, 0],
            data.psi[it, 1, :, 0], vmin = - umax, vmax = umax,
            shading = "gouraud", cmap = "seismic")
    meridional_wind = axes[0, 1].pcolormesh(data.xx[:, 0], data.zz[:, 0],
            data.psi[it, 2, :, 0], vmin = - vmax, vmax = vmax,
            shading = "gouraud", cmap = "seismic")
    vertical_wind = axes[0, 2].pcolormesh(data.xx[:, 0], data.zz[:, 0],
            data.psi[it, 3, :, 0], vmin = - wmax, vmax = wmax,
            shading = "gouraud", cmap = "seismic")
    density = axes[1, 0].pcolormesh(data.xx[:, 0], data.zz[:, 0],
            data.psi[it, 0, :, 0], vmin = - rhopmax, vmax = rhopmax,
            shading = "gouraud", cmap = "seismic")
    exner_pressure = axes[1, 1].pcolormesh(data.xx[:, 0], data.zz[:, 0],
            data.psi[it, 4, :, 0], vmin = - pipmax, vmax = pipmax,
            shading = "gouraud", cmap = "seismic")
    potential_temperature = axes[1, 2].pcolormesh(data.xx[:, 0],
            data.zz[:, 0], data.psi[it, 5, :, 0], shading = "gouraud",
            cmap = "plasma")
    axes[0, 0].set_title("Zonal wind")
    axes[0, 1].set_title("Meridional wind")
    axes[0, 2].set_title("Vertical wind")
    axes[1, 0].set_title("Density")
    axes[1, 1].set_title("Exner pressure")
    axes[1, 2].set_title("Potential temperature")
    figure.colorbar(zonal_wind, ax = axes[0, 0])
    figure.colorbar(meridional_wind, ax = axes[0, 1])
    figure.colorbar(vertical_wind, ax = axes[0, 2])
    figure.colorbar(density, ax = axes[1, 0])
    figure.colorbar(exner_pressure, ax = axes[1, 1])
    figure.colorbar(potential_temperature, ax = axes[1, 2])
    for nn in range(6):
        axes.flat[nn].set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
        axes.flat[nn].set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
    figure.savefig("../results/mountainwave_msgwam.pdf")
    figure.savefig("../results/mountainwave_msgwam.png")
else:
    wmax = numpy.max(numpy.absolute(data.psi[it, 3, :, 0]))
    figure, axes = pyplot.subplots()
    plot = axes.pcolormesh(data.xx[:, 0], data.zz[:, 0],
            data.psi[it, 3, :, 0], vmin = - wmax, vmax = wmax,
            shading = "gouraud", cmap = "seismic")
    axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
    axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
    figure.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
    figure.savefig("../results/mountainwave_msgwam.pdf")
    figure.savefig("../results/mountainwave_msgwam.png", dpi = 500)

# Make difference plot.
if (data.psi != reference.psi).all():
    if plot_all:
        umax = numpy.max(numpy.absolute(deltapsi[it, 1, :, 0]))
        vmax = numpy.max(numpy.absolute(deltapsi[it, 2, :, 0]))
        wmax = numpy.max(numpy.absolute(deltapsi[it, 3, :, 0]))
        rhopmax = numpy.max(numpy.absolute(deltapsi[it, 0, :, 0]))
        pipmax = numpy.max(numpy.absolute(deltapsi[it, 4, :, 0]))
        thetamax = numpy.max(numpy.absolute(deltapsi[it, 5, :, 0]))
        figure, axes = pyplot.subplots(2, 3, figsize = (8.0, 6.0))
        zonal_wind = axes[0, 0].pcolormesh(data.xx[:, 0], data.zz[:, 0],
                deltapsi[it, 1, :, 0], vmin = - umax, vmax = umax,
                shading = "gouraud", cmap = "seismic")
        meridional_wind = axes[0, 1].pcolormesh(data.xx[:, 0],
                data.zz[:, 0], deltapsi[it, 2, :, 0], vmin = - vmax,
                vmax = vmax, shading = "gouraud", cmap = "seismic")
        vertical_wind = axes[0, 2].pcolormesh(data.xx[:, 0],
                data.zz[:, 0], deltapsi[it, 3, :, 0], vmin = - wmax,
                vmax = wmax, shading = "gouraud", cmap = "seismic")
        density = axes[1, 0].pcolormesh(data.xx[:, 0], data.zz[:, 0],
                deltapsi[it, 0, :, 0], vmin = - rhopmax, vmax = rhopmax,
                shading = "gouraud", cmap = "seismic")
        exner_pressure = axes[1, 1].pcolormesh(data.xx[:, 0],
                data.zz[:, 0], deltapsi[it, 4, :, 0], vmin = - pipmax,
                vmax = pipmax, shading = "gouraud", cmap = "seismic")
        potential_temperature = axes[1, 2].pcolormesh(data.xx[:, 0],
                data.zz[:, 0], deltapsi[it, 5, :, 0], vmin = - thetamax,
                shading = "gouraud", cmap = "seismic")
        axes[0, 0].set_title("Zonal wind")
        axes[0, 1].set_title("Meridional wind")
        axes[0, 2].set_title("Vertical wind")
        axes[1, 0].set_title("Density")
        axes[1, 1].set_title("Exner pressure")
        axes[1, 2].set_title("Potential temperature")
        figure.colorbar(zonal_wind, ax = axes[0, 0])
        figure.colorbar(meridional_wind, ax = axes[0, 1])
        figure.colorbar(vertical_wind, ax = axes[0, 2])
        figure.colorbar(density, ax = axes[1, 0])
        figure.colorbar(exner_pressure, ax = axes[1, 1])
        figure.colorbar(potential_temperature, ax = axes[1, 2])
        for nn in range(6):
            axes.flat[nn].set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
            axes.flat[nn].set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
        figure.savefig("../results/mountainwave_msgwam_difference.pdf")
        figure.savefig("../results/mountainwave_msgwam_difference.png")
    else:
        wmax = numpy.max(numpy.absolute(deltapsi[it, 3, :, 0]))
        figure, axes = pyplot.subplots()
        plot = axes.pcolormesh(data.xx[:, 0], data.zz[:, 0],
                deltapsi[it, 3, :, 0], vmin = - wmax, vmax = wmax,
                shading = "gouraud", cmap = "seismic")
        axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
        axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
        figure.colorbar(plot, label = r"$\Delta w \,"
                r"\mathrm{\left[m \, s^{- 1}\right]}$")
        figure.savefig("../results/mountainwave_msgwam_difference.pdf")
        figure.savefig("../results/mountainwave_msgwam_difference.png",
                dpi = 500)
