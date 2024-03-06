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
  data_path = "/scratch/atmodynamics/" + user_name + "/PF/runs"
  reference_path = "/scratch/atmodynamics/" + user_name + "/PF/pinc/reference"

else:
  # Local machine
  data_path = ".."
  reference_path = data_path

print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + "/IGW/")
reference = tools.ModelOutput(reference_path + "/IGW/")

# Set plot parameters.
choice = "xz"
it = 1
ix = 32
iy = 0
iz = 14

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

# Compute difference.
deltapsi = psi - psiref

# Make plot.
maximum = numpy.max(numpy.abs(psi[5]))
figure, axes = pyplot.subplots()
plot = axes.pcolormesh(xx, yy, psi[5], vmin = - maximum, vmax = maximum, \
    shading = "gouraud", cmap = "seismic")
axes.set_xlabel(xlabel)
axes.set_ylabel(ylabel)
figure.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
figure.savefig("".join((data_path, "/results/IGW.pdf")))
figure.savefig("".join((data_path, "/results/IGW.png")), dpi = 500)

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
  figure.savefig("".join((data_path, "/results/IGW_difference.pdf")))
  figure.savefig("".join((data_path, "/results/IGW_difference.png")), dpi = 500)
