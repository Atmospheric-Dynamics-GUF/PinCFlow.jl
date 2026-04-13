# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples test/plot_gif_ags.jl

using HDF5
using CairoMakie

run = "2202_01"

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

q_obs  = Observable(data["q"][:, iy, :, 1])
qv_obs = Observable(data["qv"][:, iy, :, 1])
n_obs  = Observable(data["n"][:, iy, :, 1])
w_obs  = Observable(data["w"][:, iy, :, 1])

q_vals  = data["q"][:, iy, :, 1:tmax]
qv_vals = data["qv"][:, iy, :, 1:tmax]
w_vals  = data["w"][:, iy, :, 1:tmax]

n_vals = data["n"][:, iy, :, 1:tmax]

#println("q =", q_vals)

q_min  = minimum(skipmissing(filter(isfinite, q_vals)))
q_max  = maximum(skipmissing(filter(isfinite, q_vals)))

println("q min/max  = ", q_min, " ", q_max)

qv_min = minimum(skipmissing(filter(isfinite, qv_vals)))
qv_max = maximum(skipmissing(filter(isfinite, qv_vals)))

println("qv min/max = ", qv_min, " ", qv_max)

n_min = minimum(skipmissing(filter(isfinite, n_vals)))
n_max = maximum(skipmissing(filter(isfinite, n_vals)))

println("n min/max  = ", n_min, " ", n_max)

w_min = minimum(skipmissing(filter(isfinite, w_vals)))
w_max = maximum(skipmissing(filter(isfinite, w_vals)))

# use symmetric range for w
w_abs_max = max(abs(w_min), abs(w_max))
w_min = -w_abs_max
w_max = w_abs_max

println("w min/max  = ", w_min, " ", w_max)

fig = Figure(size = (1200, 800))

ax1 = Axis(fig[1,1], xlabel="x (km)", ylabel="z (km)", title="q")
ax2 = Axis(fig[2,1], xlabel="x (km)", ylabel="z (km)", title="qv")
ax3 = Axis(fig[1,3], xlabel="x (km)", ylabel="z (km)", title="n")
ax4 = Axis(fig[2,3], xlabel="x (km)", ylabel="z (km)", title="w")

hm1 = heatmap!(ax1, x, z, q_obs; colormap=:blues, colorrange=(q_min, q_max))
hm2 = heatmap!(ax2, x, z, qv_obs; colormap=:blues, colorrange=(qv_min, qv_max))
hm3 = heatmap!(ax3, x, z, n_obs; colormap=:blues, colorrange=(n_min, n_max))
hm4 = heatmap!(ax4, x, z, w_obs; colormap=:seismic, colorrange=(w_min, w_max))

Colorbar(fig[1,2], hm1, label="q [kg/kg]")
Colorbar(fig[2,2], hm2, label="qv [kg/kg]")
Colorbar(fig[1,4], hm3, label="n [#/kg]")
Colorbar(fig[2,4], hm4, label="w [m/s]")    
record(fig,
       output_file,
       1:tmax;
       framerate = 10) do t

    q_obs[]  = data["q"][:, iy, :, t]
    qv_obs[] = data["qv"][:, iy, :, t]
    n_obs[]  = data["n"][:, iy, :, t]
    w_obs[]  = data["w"][:, iy, :, t]
    
    ax1.title = "q at t = $t"
    ax2.title = "qv at t = $t"
    ax3.title = "n at t = $t"
    ax4.title = "w at t = $t"
end

close(data)

println("GIF saved to $output_file")