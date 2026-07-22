# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples test/fast_plot_q.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

# Base font-size (points). Change this to scale all plot text.
FS = 12

tau_vals = [0.0, 3.0e-8, 0.0003, 0.03]

data_no_sink = h5open("../ice_mountain_wave_with_smaller_sink_$(tau_vals[1]).h5", "r")
data_very_large_sink = h5open("../added_n_sink/ice_mountain_wave_tau_$(tau_vals[2])_period_1800.0.h5", "r")
#data_large_sink = h5open("../ice_mountain_wave_with_smaller_sink_$(tau_vals[3]).h5", "r")
#data_small_sink = h5open("../ice_mountain_wave_with_smaller_sink_$(tau_vals[4]).h5", "r")
#data = h5open("../u/ice_mountain_wave_u.h5", "r")
# Set grid for data1
x = data_no_sink["x"][:] .* 0.001 
y = data_no_sink["y"][:] .* 0.001
z = data_no_sink["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Set grid for data2
# x2 = data2s["x2"][:] .* 0.001 
# y2 = data2s["y2"][:] .* 0.001 
# z2 = data2s["z2"][:, :, :] .* 0.001
# x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
# y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

#iy = Int(length(data["w"][1, :, 1, 1]) / 2)
ny = size(data_no_sink["w"], 2)
iy = cld(ny, 2)
iys = 1 
println("iy = ", iy)

tidx = length(data_no_sink["w"][1, 1, 1, :])
#tidx2s = length(data2s["w"][1, 1, 1, :])
println("tidx = ", tidx)

# time index where q becomes nan somewhere
#tidx_nan_q = 0
#for t in 1:tidx
#    q_slice = data_no_sink["q"][:, :, :, t]
#    if any(isnan, q_slice)
#        tidx_nan_q = t
#        break
#    end
#end
#println("First time index with NaN in q: ", tidx_nan_q)

#fld_q_no_sink = data_no_sink["q"][:, :, :, tidx]
#fld_q_small_sink = data_small_sink["q"][:, :, :, tidx]
#fld_q_large_sink = data_large_sink["q"][:, :, :, tidx]
fld_q_very_large_sink = data_very_large_sink["q"][:, :, :, tidx]
fld_n_very_large_sink = data_very_large_sink["n"][:, :, :, tidx]
fld_q_no_sink = data_no_sink["q"][:, :, :, tidx]
#fld2s = data_small_sink["sn"][:, :, :, tidx2s]
println("size fld  ", size(fld_q_no_sink))
#println("size fld2s", size(fld2s))
println("fld  ", maximum(fld_q_no_sink), " ", minimum(fld_q_no_sink))
##println("fld2s", maximum(fld2s), " ", minimum(fld2s))

# Plot comparison: data vs datas
figure(figsize=(12, 10))
subplot(511)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_q_no_sink[:,iy,:]; cmap="Blues", label="q [kg/kg]")
cbar = colorbar(contour, label="q [kg/kg]")
cbar.ax.tick_params(labelsize=FS)
title("q no sink", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
ax = gca(); ax.tick_params(labelsize=FS)

subplot(512)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_q_very_large_sink[:,iy,:]; cmap="Blues", label="q [kg/kg]")
cbar = colorbar(contour, label="q [kg/kg]")
cbar.ax.tick_params(labelsize=FS)
title("q tau = $(tau_vals[2])", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
## ax = gca(); ax.tick_params(labelsize=FS)

subplot(513)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_q_very_large_sink[:,iy,:]; cmap="Blues", label="q [kg/kg]")
cbar = colorbar(contour, label="q [kg/kg]")
cbar.ax.tick_params(labelsize=FS)
title("q tau = $(tau_vals[2])", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
## ax = gca(); ax.tick_params(labelsize=FS)

subplot(514)
diff_q = fld_q_no_sink - fld_q_very_large_sink
contour = pcolormesh(x[:, iy, :], z[:, iy, :], diff_q[:,iy,:]; cmap="bwr", label="q difference [kg/kg]", vmin=-maximum(abs.(diff_q)), vmax=maximum(abs.(diff_q)))
cbar = colorbar(contour, label="q difference [kg/kg]")
cbar.ax.tick_params(labelsize=FS)
title("q no sink - q tau = $(tau_vals[2])", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
## ax = gca(); ax.tick_params(labelsize=FS)

subplot(515)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld_n_very_large_sink[:,iy,:]; cmap="Blues", label="n [#/kg]")
cbar = colorbar(contour, label="n [#/kg]")
cbar.ax.tick_params(labelsize=FS)
title("n tau = $(tau_vals[2])", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
## ax = gca(); ax.tick_params(labelsize=FS)

tight_layout()
savefig("../added_n_sink/output/ice_mountain_2D_with_sinks_$(tau_vals[1])_$(tau_vals[2]).png")
