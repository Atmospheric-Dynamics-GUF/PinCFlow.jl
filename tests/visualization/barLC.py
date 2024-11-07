import subprocess
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
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
data = tools.ModelOutput(data_path + 'barLC/')
reference = tools.ModelOutput(reference_path + 'barLC/')

# Set slices.
it = - 1
iz = 2
iy = slice(42, None)

print(" ".join(("Time:", str(data.variables['t'][it]), "s")))

# Set grid.
xx = data.variables['x'][:] * 0.001
yy = data.variables['y'][iy] * 0.001
dx = np.abs(xx[1] - xx[0]) * 1000.0
dy = np.abs(yy[1] - yy[0]) * 1000.0

# Loop over data and reference.
for data_set in (reference, data):

  # Set fields of interest.
  uu = data_set.groups['atmvar'].variables['u'][it, iz, iy]
  vv = data_set.groups['atmvar'].variables['v'][it, iz, iy]
  thetap = data_set.groups['atmvar'].variables['thetap'][it, iz, iy]

  # Compute divergence.
  divergence = np.zeros_like(uu)
  divergence[1:, 1:] = (uu[1:, 1:] - uu[1:, :(- 1)]) / dx + (vv[1:, 1:] - \
      vv[:(- 1), 1:]) / dy

  # Apply Fourier filter.
  urossby = np.zeros_like(uu)
  ugravity = np.zeros_like(uu)
  sigma = divergence
  sigmatilde = fft.fft2(sigma)
  kk = fft.fftfreq(sigma.shape[1], d = dx)
  ll = fft.fftfreq(sigma.shape[0], d = dy)
  urossbytilde = sigmatilde.copy()
  urossbytilde[np.abs(ll) > 1.0e-6, :] = 0.0
  urossbytilde[:, np.abs(kk) > 1.0e-6] = 0.0
  urossby = fft.ifft2(urossbytilde).real
  ugravity = sigma - urossby

  # Compute difference of relevant fields.
  if data_set == reference:
    deltaugravity = ugravity.copy()
    deltatheta = thetap.copy()
  elif data_set == data:
    deltaugravity = ugravity - deltaugravity
    deltatheta = thetap - deltatheta

# Create plot.
(levels, colormap) = style.symmetric_contours(ugravity.min(), ugravity.max())
fig, axes = plt.subplots()
plot = axes.contourf(xx, yy, ugravity, levels, cmap = colormap)
axes.contour(xx, yy, thetap, linewidths = 1.0, colors = "black")
axes.set_xlabel(r"$x \, \left[\mathrm{km}\right]$")
axes.set_ylabel(r"$y \, \left[\mathrm{km}\right]$")
fig.colorbar(plot, label = r"$\boldsymbol{\nabla}_z \cdot" r"\boldsymbol{u}" \
    r"\, \left[\mathrm{s^{- 1}}\right]$")
fig.savefig("".join((data_path, "/results/barLC.pdf")))
fig.savefig("".join((data_path, "/results/barLC.png")), dpi = 500)

# Create difference plot.
if data_path != reference_path:
  (levels, colormap) = style.symmetric_contours(deltaugravity.min(), \
      deltaugravity.max())
  fig, axes = plt.subplots()
  plot = axes.contourf(xx, yy, deltaugravity, levels, cmap = colormap)
  axes.contour(xx, yy, deltatheta, linewidths = 1.0, colors = "black")
  axes.set_xlabel(r"$x \, \left[\mathrm{km}\right]$")
  axes.set_ylabel(r"$y \, \left[\mathrm{km}\right]$")
  fig.colorbar(plot, label = r"$\Delta \boldsymbol{\nabla}_z \cdot" \
      r"\boldsymbol{u} \, \left[\mathrm{s^{- 1}}\right]$")
  fig.savefig("".join((data_path, "/results/barLC_difference.pdf")))
  fig.savefig("".join((data_path, "/results/barLC_difference.png")), \
      dpi = 500)