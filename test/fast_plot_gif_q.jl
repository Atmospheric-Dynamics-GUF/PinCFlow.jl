# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples test/fast_plot_gif_q.jl

using HDF5
using CairoMakie

run = "2101_01"

input_file = "/work/bb1097/b383844/PinCFlow/qv_forcing/results/ice_mountain_wave_$(run).h5"

data = h5open(input_file, "r")

output_file = "../qv_forcing/visualization/ice_mountain_2D_$(run)_tau.gif"

x = data["x"][:] .* 0.001
z = data["z"][1, 1, :] .* 0.001

ny = size(data["q"], 2)
iy = cld(ny, 2)
tmax = size(data["q"], 4)

println("tmax = ", tmax)

q_obs  = Observable(data["q"][:, iy, :, 1])
qv_obs = Observable(data["qv"][:, iy, :, 1])
n_obs  = Observable(data["n"][:, iy, :, 1])
dqv_obs = Observable(data["iaux3"][:, iy, :, 1])

q_vals  = data["q"][:, iy, :, :]
qv_vals = data["qv"][:, iy, :, :]
dqv_vals = data["iaux3"][:, iy, :, :]

n_vals = data["n"][:, iy, :, :]

#println("q =", q_vals)

q_min  = minimum(skipmissing(filter(isfinite, q_vals)))
q_max  = maximum(skipmissing(filter(isfinite, q_vals)))

println("q min/max  = ", q_min, " ", q_max)

qv_min = minimum(skipmissing(filter(isfinite, qv_vals)))
qv_max = maximum(skipmissing(filter(isfinite, qv_vals)))

println("qv min/max = ", qv_min, " ", qv_max)

dqv_min = minimum(skipmissing(filter(isfinite, dqv_vals)))
dqv_max = maximum(skipmissing(filter(isfinite, dqv_vals)))

println("dqv min/max = ", dqv_min, " ", dqv_max)

n_min = minimum(skipmissing(filter(isfinite, n_vals)))
n_max = maximum(skipmissing(filter(isfinite, n_vals)))

println("n min/max  = ", n_min, " ", n_max)

qv_diff = data["qv"][:, iy, :, tmax] .- data["qv"][:, iy, :, 1]

fig = Figure(size = (800, 1000))

ax1 = Axis(fig[1,1], xlabel="x (km)", ylabel="z (km)", title="q")
ax2 = Axis(fig[2,1], xlabel="x (km)", ylabel="z (km)", title="qv")
ax3 = Axis(fig[3,1], xlabel="x (km)", ylabel="z (km)", title="dqv")
ax4 = Axis(fig[4,1], xlabel="x (km)", ylabel="z (km)", title="n")
ax5 = Axis(fig[5,1], xlabel="x (km)", ylabel="z (km)", title="qv difference t=tmax - t=0")

hm1 = heatmap!(ax1, x, z, q_obs; colormap=:blues, colorrange=(q_min, q_max))
hm2 = heatmap!(ax2, x, z, qv_obs; colormap=:blues, colorrange=(qv_min, qv_max))
hm3 = heatmap!(ax3, x, z, dqv_obs; colormap=:blues, colorrange=(dqv_min, dqv_max))
hm4 = heatmap!(ax4, x, z, n_obs; colormap=:blues, colorrange=(n_min, n_max))
hm5 = heatmap!(ax5, x, z, qv_diff; colormap=:blues)

Colorbar(fig[1,2], hm1, label="q [kg/kg]")
Colorbar(fig[2,2], hm2, label="qv [kg/kg]")
Colorbar(fig[3,2], hm3, label="dqv [kg/kg/s]")
Colorbar(fig[4,2], hm4, label="n [#/kg]")
Colorbar(fig[5,2], hm5, label="qv difference [kg/kg]")

record(fig,
       output_file,
       1:tmax;
       framerate = 10) do t

    q_obs[]  = data["q"][:, iy, :, t]
    qv_obs[] = data["qv"][:, iy, :, t]
    dqv_obs[] = data["iaux3"][:, iy, :, t]
    n_obs[]  = data["n"][:, iy, :, t]
    
    ax1.title = "q at t = $t"
    ax2.title = "qv at t = $t"
    ax3.title = "dqv at t = $t"
    ax4.title = "n at t = $t"
end

close(data)

println("GIF saved to $output_file")