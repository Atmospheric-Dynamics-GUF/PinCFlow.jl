# RUN in PF/pinc/
# julia --project=. examples/visualization/fast_plot.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

#data = h5open("/home/dolaptch/PF/pinc/test/pincflow_output.h5")
#data = h5open("/home/dolaptch/PF/runs/tjl05/pincflow_output.h5")
data = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl05/pincflow_output.h5")
#data2 = h5open("/home/dolaptch/PF/runs/tjl06/pincflow_output.h5")
data2 = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5")

# Set grid for data1
x = data["x"][:] .* 0.001 .- 10
y = data["y"][:] .* 0.001 .- 10
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Set grid for data2
x2 = data2["x2"][:] .* 0.001 .- 10
y2 = data2["y2"][:] .* 0.001 .- 10
z2 = data2["z2"][:, :, :] .* 0.001
x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

iy = 1
tidx = 1 #length(data["w"][1, 1, 1, :])
tidx2 = length(data2["w"][1, 1, 1, :])

fld = data["w"][:, :, :, tidx]
fld2 = data2["sn"][:, :, :, tidx2]

println("size fld ", size(fld))
println("size fld2 ", size(fld2))
println("fld  ", maximum(fld), " ", minimum(fld))
println("fld2 ", maximum(fld2), " ", minimum(fld2))

# Plot both fields side by side
figure()
subplot(121)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
colorbar(contour)
title("Res")
xlabel("x (km)"); ylabel("z (km)")

subplot(122)
contour2 = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
colorbar(contour2)
title("Par")
xlabel("x (km)")


savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave.png")