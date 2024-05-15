import subprocess
import numpy
import numpy.fft as fft
import matplotlib.pyplot as pyplot
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
  data_path = ".."
  reference_path = data_path

print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + "/barLC/")
reference = tools.ModelOutput(reference_path + "/barLC/")

# Adust coordinate unit.
data.xx *= 0.001
data.yy *= 0.001

# Set fields of indices.
iz = int(0.1 * data.nz)
it = - 1

# Print time.
print(" ".join(("Time:", str(data.tt[it]), "s")))

# Loop over data and reference.
for data_set in (reference, data):

  # Set fields of interest.
  rho = data_set.psi[it, 0, iz, int(0.5 * data_set.ny):data_set.ny]
  uu = data_set.psi[it, 1, iz, int(0.5 * data_set.ny):data_set.ny]
  vv = data_set.psi[it, 2, iz, int(0.5 * data_set.ny):data_set.ny]
  ww = data_set.psi[it, 3, iz, int(0.5 * data_set.ny):data_set.ny]
  pi = data_set.psi[it, 4, iz, int(0.5 * data_set.ny):data_set.ny]
  theta = data_set.psi[it, 5, 0, int(0.5 * data_set.ny):data_set.ny]

  # Compute divergence.
  divergence = numpy.zeros_like(uu)
  divergence[1:, 1:] = (uu[1:, 1:] - uu[1:, :(- 1)]) / data_set.dx + (vv[1:, \
      1:] - vv[:(- 1), 1:]) / data_set.dy

  # Apply Fourier filter.
  urossby = numpy.zeros_like(uu)
  ugravity = numpy.zeros_like(uu)
  sigma = divergence
  sigmatilde = fft.fft2(sigma)
  kk = fft.fftfreq(sigma.shape[1], d = data_set.dx)
  ll = fft.fftfreq(sigma.shape[0], d = data_set.dy)
  urossbytilde = sigmatilde.copy()
  urossbytilde[numpy.abs(ll) > 1.0E-6, :] = 0.0
  urossbytilde[:, numpy.abs(kk) > 1.0E-6] = 0.0
  urossby = fft.ifft2(urossbytilde).real
  ugravity = sigma - urossby

  # Compute differences of relevant fields.
  if data_set == reference:
    deltaugravity = ugravity.copy()
    deltau = uu.copy()
    deltav = vv.copy()
    deltatheta = theta.copy()
  elif data_set == data:
    deltaugravity = ugravity - deltaugravity
    deltau = uu - deltau
    deltav = vv - deltav
    deltatheta = theta - deltatheta

# Adjust plotting area.
data.xx = data.xx[iz, int(0.5 * data.ny):data.ny]
data.yy = data.yy[iz, int(0.5 * data.ny):data.ny]

# Make plot.
peak = numpy.max(numpy.abs(ugravity))
thetamax = numpy.max(numpy.abs(theta))
figure, axes = pyplot.subplots()
plot = axes.pcolormesh(data.xx, data.yy, ugravity, vmax = peak, vmin = - peak, \
    shading = "gouraud", cmap = "seismic")
axes.quiver(data.xx[::3, ::3], data.yy[::3, ::3], uu[::3, ::3], vv[::3, ::3], \
    width = 0.01, scale = 500)
axes.contour(data.xx, data.yy, theta, linewidths = 1.0, colors = "black")
axes.set_xlabel(r"$x \, \left[\mathrm{km}\right]$")
axes.set_ylabel(r"$y \, \left[\mathrm{km}\right]$")
figure.colorbar(plot, label = r"$\boldsymbol{\nabla}_z \cdot" \
    r"\boldsymbol{u} \, \left[\mathrm{s^{- 1}}\right]$")
figure.savefig("".join((data_path, "/results/barLC.pdf")))
figure.savefig("".join((data_path, "/results/barLC.png")), dpi = 500)

# Make difference plot.
if data_path != reference_path:
  peak = numpy.max(numpy.abs(deltaugravity))
  thetamax = numpy.max(numpy.abs(deltatheta))
  figure, axes = pyplot.subplots()
  plot = axes.pcolormesh(data.xx, data.yy, deltaugravity, vmax = peak, vmin = \
      - peak, shading = "gouraud", cmap = "seismic")
  axes.quiver(data.xx[::3, ::3], data.yy[::3, ::3], deltau[::3, ::3], \
      deltav[::3, ::3], width = 0.01, scale = 500)
  axes.contour(data.xx, data.yy, deltatheta, linewidths = 1.0, colors = "black")
  axes.set_xlabel(r"$x \, \left[\mathrm{km}\right]$")
  axes.set_ylabel(r"$y \, \left[\mathrm{km}\right]$")
  figure.colorbar(plot, label = r"$\Delta \boldsymbol{\nabla}_z \cdot" \
      r"\boldsymbol{u} \, \left[\mathrm{s^{- 1}}\right]$")
  figure.savefig("".join((data_path, "/results/barLC_difference.pdf")))
  figure.savefig("".join((data_path, "/results/barLC_difference.png")), \
      dpi = 500)
