# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples test/plot_q.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

# Base font-size (points). Change this to scale all plot text.
FS = 12

data = h5open("../u/ice_mountain_wave_u.h5", "r")
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

#iy = Int(length(data["w"][1, :, 1, 1]) / 2)
ny = size(data["w"], 2)
iy = cld(ny, 2)
iys = 1 
println("iy = ", iy)

tidx = length(data["w"][1, 1, 1, :])
#tidx2s = length(data2s["w"][1, 1, 1, :])
println("tidx = ", tidx)


fld_q = data["q"][:, :, :, tidx]
fld_n = data["n"][:, :, :, tidx]
fld_u = data["u"][:, :, :, tidx]
#fld2s = data_small_sink["sn"][:, :, :, tidx2s]
println("size fld  ", size(fld_q))
#println("size fld2s", size(fld2s))
println("fld  ", maximum(fld_q), " ", minimum(fld_q))
##println("fld2s", maximum(fld2s), " ", minimum(fld2s))

# Plot comparison: data vs datas
figure(figsize=(12, 10))
subplot(311)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_q[:,iy,:]; cmap="Blues", label="q [kg/kg]")
cbar = colorbar(contour, label="q [kg/kg]")
title("q at y=$(y[iy,1,1]) km", fontsize=1.2*FS)
ylabel("z (km)", fontsize=FS)  

subplot(312)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_n[:,iy,:]; cmap="Blues", label="n [#/kg]")
cbar = colorbar(contour, label="n [#/kg]")
title("n at y=$(y[iy,1,1]) km", fontsize=1.2*FS)
ylabel("z (km)", fontsize=FS)  
xlabel("x (km)", fontsize=FS)

subplot(313)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_u[:,iy,:]; cmap="Blues", label="u [m/s]")
cbar = colorbar(contour, label="u [m/s]")
title("u at y=$(y[iy,1,1]) km", fontsize=1.2*FS)
ylabel("z (km)", fontsize=FS)  
xlabel("x (km)", fontsize=FS)

tight_layout()
savefig("test_output/ice_mountain_2D_u.png")
