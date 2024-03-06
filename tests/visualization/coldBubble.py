import subprocess
import numpy
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
  data_path = "/scratch/atmodynamics/" + user_name + "/output/2024/2024-03-05/tests"
  reference_path = "/scratch/atmodynamics/" + user_name + "/output/2024/2024-03-05/tests"

else:
  # Local machine
  data_path = ".."
  reference_path = data_path

print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + "/coldBubble/")
reference = tools.ModelOutput(reference_path + "/coldBubble/")

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
plot = axes.pcolormesh(xx, yy, psi[5], vmin = - maximum, vmax = maximum, \
    shading = "gouraud", cmap = "seismic")
axes.set_xlabel(xlabel)
axes.set_ylabel(ylabel)
figure.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
figure.savefig("".join((data_path, "/results/coldBubble.pdf")))
figure.savefig("".join((data_path, "/results/coldBubble.png")), dpi = 500)

# Make difference plot.
if data_path != reference_path:
  maximum = numpy.max(numpy.abs(deltapsi[5]))
  figure, axes = pyplot.subplots()
  plot = axes.pcolormesh(xx, yy, deltapsi[5], vmin = - maximum, vmax \
      = maximum, shading = "gouraud", cmap = "seismic")
  axes.set_xlabel(xlabel)
  axes.set_ylabel(ylabel)
  figure.colorbar(plot, label = r"$\Delta \theta' \," \
      r"\mathrm{\left[K\right]}$")
  figure.savefig("".join((data_path, "/results/coldBubble_difference.pdf")))
  figure.savefig("".join((data_path, "/results/coldBubble_difference.png")), \
      dpi = 500)
