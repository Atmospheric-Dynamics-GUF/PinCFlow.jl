import numpy
import matplotlib.pyplot as pyplot
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
data = tools.ModelOutput(data_path+"/IGW/")
reference = tools.ModelOutput(ref_path+"/IGW/")

# if tracer/ice2 includatmoed: 
#import sys
#file_index_opt_field = data_path+"/run00/"
#sys.path.append(file_index_opt_field)
#from index_opt_field import *
#print(iTr)

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
plot = axes.pcolormesh(xx, yy, psi[5], vmin = - maximum, vmax = maximum,
        shading = "gouraud", cmap = "seismic")
axes.set_xlabel(xlabel)
axes.set_ylabel(ylabel)
figure.colorbar(plot, label = r"$\theta' \, \mathrm{\left[K\right]}$")
figure.savefig("../results/IGW.pdf")
figure.savefig("../results/IGW.png", dpi = 500)

#print(data.psi.shape)
#print(data.psi[-1,1,1:3,0,1])
#print(reference.psi[-1,1,1:3,0,1])
#print((data.psi == reference.psi).all())
# Make difference plot.
if (not (data.psi == reference.psi).all() ):
    print('there is some difference !!!')
    maximum = numpy.max(numpy.abs(deltapsi[5]))
    figure, axes = pyplot.subplots()
    plot = axes.pcolormesh(xx, yy, deltapsi[5], vmin = - maximum,
            vmax = maximum, shading = "gouraud", cmap = "seismic")
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    figure.colorbar(plot, label = r"$\Delta \theta' \,"
            r"\mathrm{\left[K\right]}$")
    figure.savefig("../results/IGW_difference.pdf")
    figure.savefig("../results/IGW_difference.png", dpi = 500)
