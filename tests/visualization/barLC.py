import subprocess
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import netCDF4 as nc
import style

# Set script parameter.
make_animation = False

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
    reference_path = ".."

print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = nc.Dataset(data_path + '/pincflow_data.nc')
reference = nc.Dataset(reference_path + '/pincflow_data.nc')

x = data.variables['x'][:]
y = data.variables['y'][:] + 0.84e7

# Set fields of indices.
iz = int(0.1 * len(data.dimensions['z']))
it = - 1

# Print time.
print(" ".join(("Time:", str(data.variables['time'][it]), "s")))

# Loop over data and reference.
for data_set in (reference, data):

    # Set fields of interest.
  nx = len(data_set.dimensions['x'])
  ny = len(data_set.dimensions['y'])

  dx = np.abs(x[1] - x[0])
  dy = np.abs(y[1] - y[0])

  rho = data_set.groups['atmvar'].variables['rho'][it, iz, int(0.5 * ny):ny, :]
  uu = data_set.groups['atmvar'].variables['u'][it, iz, int(0.5 * ny):ny, :]
  vv = data_set.groups['atmvar'].variables['v'][it, iz, int(0.5 * ny):ny, :]
  ww = data_set.groups['atmvar'].variables['w'][it, iz, int(0.5 * ny):ny, :]
  pi = data_set.groups['atmvar'].variables['pi'][it, iz, int(0.5 * ny):ny, :]
  theta = data_set.groups['atmvar'].variables['theta'][it, 0, int(0.5 \
      * ny):ny, :]

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
    deltau = uu.copy()
    deltav = vv.copy()
    deltatheta = theta.copy()
  elif data_set == data:
    deltaugravity = ugravity - deltaugravity
    deltau = uu - deltau
    deltav = vv - deltav
    deltatheta = theta - deltatheta

xx = x * 0.001
yy = y[int(0.5 * ny):ny] * 0.001
xx, yy = np.meshgrid(xx, yy)

# Make plot.
peak = np.max(np.abs(ugravity))
thetamax = np.max(np.abs(theta))
fig, axes = plt.subplots()
plot = axes.pcolormesh(xx, yy, ugravity, vmax = peak, vmin = - peak, shading \
    = 'gouraud', cmap = 'seismic')
axes.quiver(xx[::3, ::3], yy[::3, ::3], uu[::3, ::3], vv[::3, ::3], width \
    = 0.01, scale = 500)
axes.contour(xx, yy, theta, linewidths = 1.0, colors = 'black')
axes.set_xlabel(r"$x \, \left[\mathrm{km}\right]$")
axes.set_ylabel(r"$y \, \left[\mathrm{km}\right]$")
fig.colorbar(plot, label = r"$\boldsymbol{\nabla}_z \cdot" r"\boldsymbol{u} \, \left[\mathrm{s^{- 1}}\right]$")
fig.savefig("".join((data_path, "/results/barLC.pdf")))
fig.savefig("".join((data_path, "/results/barLC.png")), dpi = 500)

if data_path != reference_path:
  peak = np.max(np.abs(deltaugravity))
  thetamax = np.max(np.abs(deltatheta))
  fig, axes = plt.subplots()
  plot = axes.pcolormesh(xx, yy, deltaugravity, vmax = peak, vmin = - peak, \
      shading = 'gouraud', cmap = 'seismic')
  axes.quiver(xx[::3, ::3], yy[::3, ::3], deltau[::3, ::3], deltav[::3, ::3], \
      width = 0.01, scale = 500)
  axes.contour(xx, yy, deltatheta, linewidths = 1.0, colors = 'black')
  axes.set_xlabel(r"$x \, \left[\mathrm{km}\right]$")
  axes.set_ylabel(r"$y \, \left[\mathrm{km}\right]$")
  fig.colorbar(plot, label = r"$\Delta \boldsymbol{\nabla}_z \cdot" r"\boldsy" \
      r"mbol{u} \, \left[\mathrm{s^{- 1}}\right]$")
  fig.savefig("".join((data_path, "/results/barLC_difference.pdf")))
  fig.savefig("".join((data_path, "/results/barLC_difference.png")), \
      dpi = 500)