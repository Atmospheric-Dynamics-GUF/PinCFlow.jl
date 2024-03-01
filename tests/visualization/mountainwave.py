import subprocess
import numpy
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import tools
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
  reference_path = data_path

print("data_path =", data_path)
print("reference_path =", reference_path)

# Import data.
data = tools.ModelOutput(data_path + "/mountainwave/")
data.transform()
reference = tools.ModelOutput(reference_path + "/mountainwave/")
reference.transform()

# Shift grid.
data.xx -= 0.5 * data.lx
data.yy -= 0.5 * data.ly

# Adjust coordinate unit.
data.xx *= 0.001
data.yy *= 0.001
data.zz *= 0.001
data.zc *= 0.001

# Adjust output time unit.
data.tt = numpy.round(data.tt / 60.0 / 60.0, 1)

# Print last output time.
print(" ".join(("Time:", str(data.tt[- 1]), "h")))

# Set vertical velocities.
ww = data.psi[:, 3, :, 0]
wref = reference.psi[:, 3, :, 0]

# Compute difference.
deltaw = ww - wref

# Compute maximum and contour levels.
peak = numpy.max(numpy.abs(ww[- 1]))
levels = numpy.arange(- peak, peak, 0.25)

# Make plot.
figure, axes = pyplot.subplots()
plot = axes.pcolormesh(data.xx[:, 0], data.zc[- 1, :, 0], ww[- 1], cmap \
    = "seismic", vmin = - peak, vmax = peak, shading = "gouraud")
axes.contour(data.xx[:, 0], data.zc[- 1, :, 0], ww[- 1], levels = levels, \
    colors = "black", linewidths = 1.0)
axes.plot(data.xx[0, 0], 0.001 * data.hh[- 1, 0], color = "black", linewidth \
    = 1.0)
axes.set_xlim(- 10.0, 10.0)
axes.set_ylim(0.0, 10.0)
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
figure.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
figure.savefig("".join((data_path, "/results/mountainwave.pdf")))
figure.savefig("".join((data_path, "/results/mountainwave.png")), dpi = 500)

# Make animation.
if make_animation:

  def animate(it):
    it = int(0.01 * data.nt * it)
    axes.clear()
    axes.pcolormesh(data.xx[:, 0], data.zc[it, :, 0], ww[it], cmap \
        = "seismic", vmin = - peak, vmax = peak, shading = "gouraud")
    axes.contour(data.xx[:, 0], data.zc[it, :, 0], ww[it], levels = levels, \
        colors = "black", linewidths = 1.0)
    axes.plot(data.xx[0, 0], 0.001 * data.hh[it, 0], color = "black", \
        linewidth = 1.0)
    axes.set_title("".join((r"Time: $", str(data.tt[it]), r"\, \mathrm{h}$")))
    axes.set_xlim(- 10.0, 10.0)
    axes.set_ylim(0.0, 10.0)
    axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
    axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
    return axes

  figure, axes = pyplot.subplots()
  plot = axes.pcolormesh(data.xx[:, 0], data.zc[- 1, :, 0], ww[- 1], cmap \
      = "seismic", vmin = - peak, vmax = peak, shading = "gouraud")
  figure.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
  animation_output = animation.FuncAnimation(figure, animate)
  writer = animation.FFMpegWriter(fps = 10)
  animation_output.save("".join((data_path, "/results/mountainwave.avi")), \
      writer = writer)

# Make difference plot.
if data_path != reference_path:
  peak = numpy.max(numpy.abs(deltaw[- 1]))
  figure, axes = pyplot.subplots()
  plot = axes.pcolormesh(data.xx[:, 0], data.zc[- 1, :, 0], deltaw[- 1], cmap \
      = "seismic", vmin = - peak, vmax = peak, shading = "gouraud")
  axes.contour(data.xx[:, 0], data.zc[- 1, :, 0], deltaw[- 1], colors \
      = "black", linewidths = 1.0)
  axes.plot(data.xx[0, 0], 0.001 * data.hh[- 1, 0], color = "black", linewidth \
      = 1.0)
  axes.set_xlim(- 10.0, 10.0)
  axes.set_ylim(0.0, 10.0)
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
  figure.colorbar(plot, label = r"$w \, \mathrm{\left[m \, s^{- 1}\right]}$")
  figure.savefig("".join((data_path, "/results/mountainwave.pdf")))
  figure.savefig("".join((data_path, "/results/mountainwave.png")), dpi = 500)
