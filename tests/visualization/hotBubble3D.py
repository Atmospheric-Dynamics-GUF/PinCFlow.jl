import numpy
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import tools
import subprocess

name_host=subprocess.getoutput('hostname')
name_user=subprocess.getoutput('whoami')

if ('levante' in name_host):
    #Levante cluster
    data_path = "/scratch/b/"+name_user+"/PF/runs"
    ref_path = "/scratch/b/"+name_user+"/PF/pinc/reference"

elif ('login' in name_host):
    #Goethe cluster
    data_path = "/scratch/atmodynamics/"+name_user+"/PF/runs"
    ref_path = "/scratch/atmodynamics/"+name_user+"/PF/pinc/reference"

else :
    print('unknown host !')
    print('set variables data_path, ref_path')
    exit()

print("data_path=", data_path)
print("ref_path=", ref_path)    

# Import data.
data = tools.ModelOutput(data_path+"/hotBubble3D/")
reference = tools.ModelOutput(ref_path+"/hotBubble3D/")

# Adjust coordinate unit.
data.xx *= 0.001
data.zz *= 0.001

# Set y index.
it = int(0.4 * data.nt)
iy = int(0.5 * data.ny)

# Print time.
print(" ".join(("Time:", str(data.tt[it]), "s")))

# Set potential temperatures.
theta = data.psi[:, 5]
thetaref = reference.psi[:, 5]

# Compute difference.
deltatheta = theta - thetaref

# Subtract offset.
theta = theta - 300.0

# Define animation function.
def animate(nn):
    axes.clear()
    axes.pcolormesh(data.xx[:, iy], data.zz[:, iy], theta[nn, :, iy],
            cmap = "seismic", vmin = - peak, vmax = peak, shading = "gouraud")
    axes.set_title("".join((r"Time: $" + str(int(data.tt[nn])), r"\,"
            r"\mathrm{s}$")))
    axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
    axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
    return axes

# Make plot.
figure, axes = pyplot.subplots()
peak = numpy.max(numpy.absolute(theta[it, :, iy]))
plot = axes.pcolormesh(data.xx[:, iy], data.zz[:, iy],
        theta[it, :, iy], cmap = "seismic", vmin = - peak,
        vmax = peak, shading = "gouraud")
axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
figure.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
figure.savefig("../results/hotBubble3D.pdf")
figure.savefig("../results/hotBubble3D.png", dpi = 500)

# Make animation.
if make_animation:
    animation_output = animation.FuncAnimation(figure, animate,
            frames = data.nt)
    writer = animation.FFMpegWriter(fps = 10)
    animation_output.save("../results/hotBubble3D.avi", writer = writer)

# Make difference plot.
if (not (data.psi != reference.psi).all()):
    figure, axes = pyplot.subplots()
    peak = numpy.max(numpy.absolute(deltatheta[it, :, iy]))
    plot = axes.pcolormesh(data.xx[:, iy], data.zz[:, iy],
            deltatheta[it, :, iy], cmap = "seismic",vmin = - peak,
                    vmax = peak, shading = "gouraud")
    axes.set_xlabel(r"$x \, \mathrm{\left[km\right]}$")
    axes.set_ylabel(r"$z \, \mathrm{\left[km\right]}$")
    figure.colorbar(plot, label = r"$\Delta \theta' \,"
            r"\mathrm{\left[K\right]}$")
    figure.savefig("../results/hotBubble3D_difference.pdf")
    figure.savefig("../results/hotBubble3D_difference.png", dpi = 500)
