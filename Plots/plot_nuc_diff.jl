# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples test/plot_nuc_diff.jl

using HDF5
using CairoMakie

run = "0904_04"

directory = "adv"

input_file = "/work/bb1097/b383844/PinCFlow/$directory/results/ice_mountain_wave_$(run).h5"

data = h5open(input_file, "r")

output_file = "../$directory/visualization/ice_mountain_2D_$(run)_tau.mp4"
#output_file = "/work/bb1097/b383844/PinCFlow/$directory/ice_mountain_2D_$(run)_tau.gif"

x = data["x"][:] .* 0.001
z = data["z"][1, 1, :] .* 0.001

ny = size(data["q"], 2)
iy = cld(ny, 2)
tmax = size(data["q"], 4)

tmax = tmax

println("tmax = ", tmax)

n_obs = Observable(data["n"][:, iy, :, 1])
n_vals = data["n"][:, iy, :, 1:tmax]
n_Nuc_vals = data["nNuc"][:, iy, :, 1:tmax]
n_diff_obs = Observable(n_Nuc_vals[:, :, 2] - n_Nuc_vals[:, :, 1])
n_diff_vals = n_Nuc_vals[:,:,2:tmax] - n_Nuc_vals[:,:,1:tmax-1]

n_min = minimum(skipmissing(filter(isfinite, n_vals)))
n_max = maximum(skipmissing(filter(isfinite, n_vals)))

n_diff_min = minimum(skipmissing(filter(isfinite, n_diff_vals)))
n_diff_max = maximum(skipmissing(filter(isfinite, n_diff_vals)))

println("n range: ", n_min, " to ", n_max)
println("n diff range: ", n_diff_min, " to ", n_diff_max)

fig = Figure(size = (800, 600))

ax1 = Axis(fig[1,1], xlabel="x (km)", ylabel="z (km)", title="n")
ax2 = Axis(fig[2,1], xlabel="x (km)", ylabel="z (km)", title="n diff")

hm1 = heatmap!(ax1, x, z, n_obs; colormap=:blues, colorrange=(n_min, n_max))
hm2 = heatmap!(ax2, x, z, n_diff_obs; colormap=:blues, colorrange=(n_diff_min, n_diff_max))

Colorbar(fig[1,2], hm1, label="n [#/kg]")
Colorbar(fig[2,2], hm2, label="n diff [#/kg]")

record(fig, output_file, 2:tmax; framerate=10) do t
    n_obs[] = data["n"][:, iy, :, t]
    n_diff_obs[] = n_Nuc_vals[:, :, t] - n_Nuc_vals[:, :, t-1]

    ax1.title = "n at t = $t"
    ax2.title = "n diff at t = $t"
end

close(data)

println("GIF saved to $output_file")