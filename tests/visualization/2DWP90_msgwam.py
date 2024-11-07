import subprocess
import numpy as np
import matplotlib.pyplot as plt
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
data = tools.ModelOutput(data_path + '2DWP90_msgwam/')
reference = tools.ModelOutput(reference_path + '2DWP90_msgwam/')

# Set slices.
it = - 1
iz = slice(15, 45)

print(" ".join(("Time:", str(data.variables['t'][it]), "s")))

# Set grid.
xx = data.variables['x'][:] * 0.001
zz = data.variables['z'][iz] * 0.001

# Get large-scale tracer distribution change (tLS)
# and zonal, meridional, and vertical wind (uLS, vLS, wLS)
tLS = data.groups['atmvar'].variables['tmrd'][it, iz, 0]
uLS = data.groups['atmvar'].variables['u'][it, iz, 0]
vLS = data.groups['atmvar'].variables['v'][it, iz, 0]
wLS = data.groups['atmvar'].variables['w'][it, iz, 0]

tLSref = reference.groups['atmvar'].variables['tmrd'][it, iz, 0]
uLSref = reference.groups['atmvar'].variables['u'][it, iz, 0]
vLSref = reference.groups['atmvar'].variables['v'][it, iz, 0]
wLSref = reference.groups['atmvar'].variables['w'][it, iz, 0]

deltat = tLS - tLSref
deltau = uLS - uLSref
deltav = vLS - vLSref
deltaw = wLS - wLSref

# Create plot.
(fig, axes) = plt.subplots(nrows = 2, ncols = 2, figsize = (8, 6))

for ax in axes.flatten():
  ax.set_xlabel(r'$x \left[\mathrm{km}\right]$')
  ax.set_ylabel(r'$z \left[\mathrm{km}\right]$')

axes[0, 0].set_title(r'Large-scale zonal wind $u \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
(levels, colormap) = style.symmetric_contours(uLS.min(), uLS.max())
pplot00 = axes[0, 0].contourf(xx, zz, uLS, levels, cmap = colormap)
fig.colorbar(pplot00, ax = axes[0, 0])

axes[0, 1].set_title(r'Large-scale meridional wind $v \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
(levels, colormap) = style.symmetric_contours(vLS.min(), vLS.max())
pplot01 = axes[0, 1].contourf(xx, zz, vLS, levels, cmap = colormap)
fig.colorbar(pplot01, ax = axes[0, 1])

axes[1, 0].set_title(r'Large-scale vertical wind $w \, \left[\mathrm{m \,' \
    r's^{- 1}}\right]$')
(levels, colormap) = style.symmetric_contours(wLS.min(), wLS.max())
pplot10 = axes[1, 0].contourf(xx, zz, wLS, levels, cmap = colormap)
fig.colorbar(pplot10, ax = axes[1, 0])

axes[1, 1].set_title('Large-scale tracer distribution')
(levels, colormap) = style.symmetric_contours(tLS.min(), tLS.max())
pplot11 = axes[1, 1].contourf(xx, zz, tLS, levels, cmap = colormap)
fig.colorbar(pplot11, ax = axes[1, 1])

fig.savefig("".join((data_path, "/results/2DWP90_msgwam.pdf")))
fig.savefig("".join((data_path, "/results/2DWP90_msgwam.png")), dpi = 500)

# Create difference plot.
if data_path != reference_path:
  fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (8, 6))

  for ax in axes.flatten():
    ax.set_xlabel(r'$x \left[\mathrm{km}\right]$')
    ax.set_ylabel(r'$z \left[\mathrm{km}\right]$')

  axes[0, 0].set_title(r'Large-scale zonal wind $u \, \left[\mathrm{m \,' \
      r's^{- 1}}\right]$')
  (levels, colormap) = style.symmetric_contours((deltau).min(), (deltau).max())
  pplot00 = axes[0, 0].contourf(xx, zz, deltau, levels, cmap = colormap)
  fig.colorbar(pplot00, ax = axes[0, 0])

  axes[0, 1].set_title(r'Large-scale meridional wind $v \, \left[\mathrm{m \,' \
      r's^{- 1}}\right]$')
  (levels, colormap) = style.symmetric_contours((deltav).min(), (deltav).max())
  pplot01 = axes[0, 1].contourf(xx, zz, deltav, levels, cmap = colormap)
  fig.colorbar(pplot01, ax = axes[0, 1])

  axes[1, 0].set_title(r'Large-scale vertical wind $w \, \left[\mathrm{m \,' \
      r's^{- 1}}\right]$')
  (levels, colormap) = style.symmetric_contours((deltaw).min(), (deltaw).max())
  pplot10 = axes[1, 0].contourf(xx, zz, deltaw, levels, cmap = colormap)
  fig.colorbar(pplot10, ax = axes[1, 0])

  axes[1, 1].set_title('Large-scale tracer distribution')
  (levels, colormap) = style.symmetric_contours((deltat).min(), (deltat).max())
  pplot11 = axes[1, 1].contourf(xx, zz, deltat, levels, cmap = colormap)
  fig.colorbar(pplot11, ax = axes[1, 1])

  fig.savefig("".join((data_path, "/results/2DWP90_msgwam_difference.pdf")))
  fig.savefig("".join((data_path, "/results/2DWP90_msgwam_difference.png")), \
      dpi = 500)