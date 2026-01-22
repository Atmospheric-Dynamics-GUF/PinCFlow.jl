# RUN in PF/pinc/
# julia --project=examples test/fast_plot_run_levante.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

# Base font-size (points). Change this to scale all plot text.
FS = 12

data_negative = h5open("ice_mountain_wave_negative_u.h5")
data_positive = h5open("ice_mountain_wave_positive_u.h5")
#data = h5open("/home/dolaptch/PF/pinc/exp/pincflow_output.h5")
#data = h5open("/home/dolaptch/PF/runs/tjl05/pincflow_output.h5")
#data = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl05/pincflow_output.h5")

#data2s = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5")
#data2s = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5")

# Set grid for data1
x = data_negative["x"][:] .* 0.001 
y = data_negative["y"][:] .* 0.001
z = data_negative["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Set grid for data2
# x2 = data2s["x2"][:] .* 0.001 
# y2 = data2s["y2"][:] .* 0.001 
# z2 = data2s["z2"][:, :, :] .* 0.001
# x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
# y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

#iy = Int(length(data["w"][1, :, 1, 1]) / 2)
ny = size(data_negative["w"], 2)
iy = cld(ny, 2)
iys = 1 
println("iy = ", iy)

tidx = length(data_negative["w"][1, 1, 1, :])
#tidx2s = length(data2s["w"][1, 1, 1, :])
println("tidx = ", tidx)


fld_negative = data_negative["q"][:, :, :, tidx]
fld_positive = data_positive["q"][:, :, :, tidx]
#fld2s = data2s["sn"][:, :, :, tidx2s]

println("size fld_negative  ", size(fld_negative))
#println("size fld2s", size(fld2s))
println("fld-  ", maximum(fld_negative), " ", minimum(fld_negative))
println("fld+ ", maximum(fld_positive), " ", minimum(fld_positive))
##println("fld2s", maximum(fld2s), " ", minimum(fld2s))

# Plot comparison: data vs datas
figure(figsize=(12, 10))
subplot(211)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_negative[:,iy,:]; cmap="Blues", label="q [kg/kg]")
cbar = colorbar(contour, label="q [kg/kg]")
cbar.ax.tick_params(labelsize=FS)
title("data: Res", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
ax = gca(); ax.tick_params(labelsize=FS)

subplot(212)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_positive[:,iy,:]; cmap="Blues", label="q [kg/kg]")
cbar = colorbar(contour, label="q [kg/kg]")
cbar.ax.tick_params(labelsize=FS)
title("data: Res", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
## ax = gca(); ax.tick_params(labelsize=FS)

tight_layout()
savefig("test_output/ice_mountain_2D_ensemble_q_with_sink.png")
