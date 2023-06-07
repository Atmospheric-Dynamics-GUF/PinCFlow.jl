import numpy
import numpy.fft as fft
import matplotlib.pyplot as pyplot
import tools
import style

data_path = "/scratch/atmodynamics/dolaptchiev/PF/pinc/tests"
ref_path = "/scratch/atmodynamics/dolaptchiev/PF/pinc/tests"

# Import data.
data = tools.ModelOutput(data_path+"/barLC/")
reference = tools.ModelOutput(ref_path+"/barLC/")

# Adust coordinate unit.
data.xx *= 0.001
data.yy *= 0.001

# Set fields of indices.
iz = int(0.1 * data.nz)
it = - 1

# Print time.
print(" ".join(("Time:", str(data.tt[it]), "s")))

# Loop over data and reference.
for set in (reference, data):

    # Set fields of interest.
    rho = set.psi[it, 0, iz, int(0.5 * set.ny):set.ny]
    uu = set.psi[it, 1, iz, int(0.5 * set.ny):set.ny]
    vv = set.psi[it, 2, iz, int(0.5 * set.ny):set.ny]
    ww = set.psi[it, 3, iz, int(0.5 * set.ny):set.ny]
    pi = set.psi[it, 4, iz, int(0.5 * data.ny):set.ny]
    rhop = set.psi[it, 5, 0, int(0.5 * set.ny):set.ny]

    # Compute divergence.
    divergence = numpy.zeros_like(uu)
    divergence[1:, 1:] = (uu[1:, 1:] - uu[1:, :(- 1)]) / set.dx + (vv[1:, 1:]
            - vv[:(- 1), 1:]) / set.dy

    # Apply Fourier filter.
    urossby = numpy.zeros_like(uu)
    ugravity = numpy.zeros_like(uu)
    sigma = divergence
    sigmatilde = fft.fft2(sigma)
    kk = fft.fftfreq(sigma.shape[1], d = set.dx)
    ll = fft.fftfreq(sigma.shape[0], d = set.dy)
    urossbytilde = sigmatilde.copy()
    urossbytilde[numpy.abs(ll) > 1.0E-6, :] = 0.0
    urossbytilde[:, numpy.abs(kk) > 1.0E-6] = 0.0
    urossby = fft.ifft2(urossbytilde).real
    ugravity = sigma - urossby

    # Compute differences of relevant fields.
    if set == reference:
        deltaugravity = ugravity.copy()
        deltau = uu.copy()
        deltav = vv.copy()
        deltarhop = rhop.copy()
    elif set == data:
        deltaugravity = ugravity - deltaugravity
        deltau = uu - deltau
        deltav = vv - deltav
        deltarhop = rhop - deltarhop

# Adjust plotting area.
data.xx = data.xx[iz, int(0.5 * data.ny):data.ny]
data.yy = data.yy[iz, int(0.5 * data.ny):data.ny]

print('Exit before plotting')

# Make plot.
peak = numpy.max(numpy.abs(ugravity))
rhopmax = numpy.max(numpy.abs(rhop))
figure, axes = pyplot.subplots()
plot = axes.pcolormesh(data.xx, data.yy, ugravity, vmax = peak, vmin = - peak,
        shading = "gouraud", cmap = "seismic")
axes.quiver(data.xx[::3, ::3], data.yy[::3, ::3], uu[::3, ::3],
        vv[::3, ::3], width = 0.01, scale = 500)
axes.contour(data.xx, data.yy, rhop,
        linewidths = 1.0, colors = "black")
print('Do plotting')

#*axes.set_xlabel(r"$x \, \left[\mathrm{km}\right]$")
#*axes.set_ylabel(r"$y \, \left[\mathrm{km}\right]$")
#*figure.colorbar(plot, label = r"$\boldsymbol{\nabla}_z \cdot"
#*        r"\boldsymbol{u} \, \left[\mathrm{s^{- 1}}\right]$")
figure.savefig("../results/barLC.png", dpi = 500)
figure.savefig("../results/barLC.pdf")

print('Plotting finished')

# Make difference plot.
if (data.psi != reference.psi).all():
    peak = numpy.max(numpy.abs(deltaugravity))
    rhopmax = numpy.max(numpy.abs(deltarhop))
    figure, axes = pyplot.subplots()
    plot = axes.pcolormesh(data.xx, data.yy, deltaugravity, vmax = peak,
            vmin = - peak, shading = "gouraud", cmap = "seismic")
    axes.quiver(data.xx[::3, ::3], data.yy[::3, ::3], deltau[::3, ::3],
            deltav[::3, ::3], width = 0.01, scale = 500)
    axes.contour(data.xx, data.yy, deltarhop,
            linewidths = 1.0, colors = "black")
#*    axes.set_xlabel(r"$x \, \left[\mathrm{km}\right]$")
#*    axes.set_ylabel(r"$y \, \left[\mathrm{km}\right]$")
#*    figure.colorbar(plot, label = r"$\Delta \boldsymbol{\nabla}_z \cdot"
#*            r"\boldsymbol{u} \, \left[\mathrm{s^{- 1}}\right]$")
#*    figure.savefig("../results/barLC_difference.pdf")
#*    figure.savefig("../results/barLC_difference.png", dpi = 500)
