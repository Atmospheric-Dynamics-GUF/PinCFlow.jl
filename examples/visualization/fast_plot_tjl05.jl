# RUN in PF/pinc/
# julia --project=. examples/visualization/fast_plot_compare.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics
using NCDatasets

data = h5open("/home/dolaptch/PF/runs/tjl05/pincflow_output.h5")  #LES
data2 = h5open("/home/dolaptch/PF/runs/tjl06/pincflow_output.h5") # RT

# Set grid.
x = data["x"][:] .* 0.001 
y = data["y"][:] .* 0.001 
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

x2 = data2["x2"][:] .* 0.001 
y2 = data2["y2"][:] .* 0.001 
z2 = data2["z2"][:, :, :] .* 0.001

x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

iy = 1
tidx = length(data["w"][1, 1, 1, :])  # for a 4D array, 4th dimension is often "time"

#fld = data["rhobar"][:, :, :]
#fld = data["thetabar"][:, :, :]
fld = data["w"][:, :, :, tidx]
#fld = data["n"][:, :, :, tidx]
fld2 = data2["sn"][:, :, :, tidx]
#fld2 = data2["wwp"][:, :, :, tidx] #+ data["w"][:, :, :, tidx]
#fld = data["thetap"][:, :, :, tidx]
#fld = data["rhop"][:, :, :, tidx]
#fld = data["n"][:, :, :, tidx]
#fld2 = data2["iaux1"][:, :, :, tidx]
#fld = data["n2"][:, :, :]
#fld0 = data["pip"][:, :, :, 1]
#fld = data["pip"][:, :, :, tidx]
#fld = data["thetabar"][:, :, :]
println("size fld ", size(fld))
println("size fld2 ", size(fld2))
println("fld  ", maximum(fld), " ", minimum(fld))
println("fld2 ", maximum(fld2), " ", minimum(fld2))

subplot(131)
subplots_adjust(wspace=1.6) 

#yaxis_limits = (7, 9.5)
yaxis_limits = (7.5, 10.5)
ylim(yaxis_limits)
#contour = contourf( fld[:,1,:]')
#contour = contourf(x2[:, iy, :], z2[:, iy, :], fld[:,iy,:]; cmap="RdBu")
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
colorbar(contour)
title(".jl")

subplot(132)
contour = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
colorbar(contour)
ylim(yaxis_limits)
title(".f90")

subplot(133)
#contour2 = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]-fld2[:,iy,:]; cmap="Blues")
#contour2 = contourf(x2[:, iy, :], z2[:, iy, :], fld[:,iy,:]-fld2[:,iy,:]; cmap="RdBu")
#colorbar(contour2)
#ylim(yaxis_limits)
title("Difference")
savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave.png")

clf()
iz = 90
plot(fld[:, iy, iz])
plot(fld2[:, iy, iz])

savefig("/home/dolaptch/PF/pinc/examples/visualization/test.png")