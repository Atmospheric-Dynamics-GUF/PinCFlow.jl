# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples Plots/test_plot.jl

""" Create test plot with n, q and w for the last time step of ice_mountain_wave_test.jl simulation. 
    The maximum values of n, q and w are displayed for comparison."""

using HDF5
using CairoMakie
using Printf

input_file = "/home/b/b383844/PinCFlow/sedimentation/results/ice_mountain_wave_test.h5"

data = h5open(input_file, "r")

output_file = "test_output/ice_mountain_wave_test.png"

x = data["x"][:] .* 0.001
z = data["z"][1, 1, :] .* 0.001

# plot n, q and w for the last time step
tmax = size(data["q"], 4)
n = data["n"][:, 1, :, tmax]
q = data["q"][:, 1, :, tmax]
w = data["w"][:, 1, :, tmax]

fig = Figure(size = (800, 1200))

ax1 = Axis(fig[1, 1], title = "Ice number concentration (n)", xlabel = "x (km)", ylabel = "z (km)")
hm1 = heatmap!(ax1, x, z, n, colormap = :blues)
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[2, 1], title = "Ice mass mixing ratio (q)", xlabel = "x (km)", ylabel = "z (km)")
hm2 = heatmap!(ax2, x, z, q, colormap = :blues)
Colorbar(fig[2, 2], hm2)

ax3 = Axis(fig[3, 1], title = "Vertical velocity (w)", xlabel = "x (km)", ylabel = "z (km)")
hm3 = heatmap!(ax3, x, z, w, colormap = :coolwarm)
Colorbar(fig[3, 2], hm3)

max_n = maximum(n)
max_q = maximum(q)
max_w = maximum(w)

text!(ax1, 0.02, 0.90, text = "Max: $(@sprintf("%.4e", max_n))", space = :relative, color = :black, fontsize = 14)
text!(ax2, 0.02, 0.90, text = "Max: $(@sprintf("%.4e", max_q))", space = :relative, color = :black, fontsize = 14)
text!(ax3, 0.02, 0.90, text = "Max: $(@sprintf("%.4e", max_w))", space = :relative, color = :black, fontsize = 14)

save(output_file, fig)

close(data)

println("Figure saved to $output_file")