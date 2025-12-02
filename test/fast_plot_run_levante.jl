# RUN in PF/pinc/
# julia --project=. examples/visualization/fast_plot.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

# Base font-size (points). Change this to scale all plot text.
FS = 12

data = h5open("mountain_wave.h5")
#data = h5open("/home/dolaptch/PF/pinc/exp/pincflow_output.h5")
#data = h5open("/home/dolaptch/PF/runs/tjl05/pincflow_output.h5")
#data = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl05/pincflow_output.h5")

#data2s = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5")
#data2s = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5")

# Set grid for data1
x = data["x"][:] .* 0.001 
y = data["y"][:] .* 0.001
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Set grid for data2
# x2 = data2s["x2"][:] .* 0.001 
# y2 = data2s["y2"][:] .* 0.001 
# z2 = data2s["z2"][:, :, :] .* 0.001
# x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
# y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

iy = Int(length(data["w"][1, :, 1, 1]) / 2)
iys = 1 
println("iy = ", iy)

tidx = length(data["w"][1, 1, 1, :])
#tidx2s = length(data2s["w"][1, 1, 1, :])


fld = data["iaux1"][:, :, :, tidx]
#fld2s = data2s["sn"][:, :, :, tidx2s]

println("size fld  ", size(fld))
#println("size fld2s", size(fld2s))
println("fld  ", maximum(fld), " ", minimum(fld))
##println("fld2s", maximum(fld2s), " ", minimum(fld2s))

# Plot comparison: data vs datas
figure(figsize=(12, 10))
subplot(211)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
cbar = colorbar(contour)
cbar.ax.tick_params(labelsize=FS)
title("data: Res", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
ax = gca(); ax.tick_params(labelsize=FS)

## subplot(212)
## contour2s = pcolormesh(x2[:, iys, :], z2[:, iys, :], fld2s[:,iys,:]; cmap="Blues")
## cbar2 = colorbar(contour2s)
## cbar2.ax.tick_params(labelsize=FS)
## title("data2s: Par (saved)", fontsize=1.5*FS)
## xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
## ax = gca(); ax.tick_params(labelsize=FS)

tight_layout()
savefig("~/output/mountain_wave.png")
