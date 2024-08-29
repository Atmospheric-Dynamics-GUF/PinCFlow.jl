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
data = tools.ModelOutput(data_path + "/hotBubble3D/")
reference = tools.ModelOutput(reference_path + "/hotBubble3D/")

# Adjust coordinate unit.
data.xx *= 0.001
data.zz *= 0.001

# Set y index.
it = int(0.4 * data.nt)
iy = int(0.5 * data.ny)

# Print time.
print(" ".join(("Time:", str(data.tt[it]), "s")))

# Set potential temperatures.
theta = data.psi["theta"]
thetaref = reference.psi["theta"]

# Compute difference.
deltatheta = theta - thetaref

# Subtract offset.
theta = theta - 300.0

# Define animation function.
def animate(nn):
  axes.clear()
  axes.pcolormesh(data.xx[:, iy], data.zz[:, iy], theta[nn, :, iy], cmap \
      = "seismic", vmin = - peak, vmax = peak, shading = "gouraud")
  axes.set_title("".join((r"Time: $" + str(int(data.tt[nn])), r"\," \
      r"\mathrm{s}$")))
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
  return axes

# Make plot.
figure, axes = pyplot.subplots()
peak = numpy.max(numpy.absolute(theta[it, :, iy]))
plot = axes.pcolormesh(data.xx[:, iy], data.zz[:, iy], theta[it, :, iy], cmap \
    = "seismic", vmin = - peak, vmax = peak, shading = "gouraud")
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
figure.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
figure.savefig("".join((data_path, "/results/hotBubble3D.pdf")))
figure.savefig("".join((data_path, "/results/hotBubble3D.png")), dpi = 500)

# Make animation.
if make_animation:
  animation_output = animation.FuncAnimation(figure, animate, frames = data.nt)
  writer = animation.FFMpegWriter(fps = 10)
  animation_output.save("".join((data_path, "/results/hotBubble3D.avi")), \
      writer = writer)

# Make difference plot.
if data_path != reference_path:
  figure, axes = pyplot.subplots()
  peak = numpy.max(numpy.absolute(deltatheta[it, :, iy]))
  plot = axes.pcolormesh(data.xx[:, iy], data.zz[:, iy], deltatheta[it, :, \
      iy], cmap = "seismic", vmin = - peak, vmax = peak, shading = "gouraud")
  axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
  axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
  figure.colorbar(plot, label = r"$\Delta \theta' \," \
      r"\mathrm{\left[K\right]}$")
  figure.savefig("".join((data_path, "/results/hotBubble3D_difference.pdf")))
  figure.savefig("".join((data_path, "/results/hotBubble3D_difference.png")), \
      dpi = 500)
