# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples test/fast_plot_ice.jl

using HDF5
using CairoMakie
using Statistics

run = "2301_14"

directory = "sedimentation"

tau_star = 3.0e-11

input_file = "/work/bb1097/b383844/PinCFlow/$directory/results/ice_mountain_wave_$(run).h5"
data = h5open(input_file, "r")

output_file = "../$directory/visualization/ice_mountain_2D_$(run)_tau.gif"

x = data["x"][:] .* 0.001
z = data["z"][1, 1, :] .* 0.001

ny = size(data["q"], 2)
iy = cld(ny, 2)
tmax = size(data["q"], 4)

println("tmax = ", tmax)

tref = 29.8 # s
mref = 7.8e11
qref = 1
nref = 1.0 / mref

q_obs = Observable(data["q"][:, iy, :, 1])
n_obs  = Observable(data["n"][:, iy, :, 1])
qv_obs = Observable(data["qv"][:, iy, :, 1])
w_obs  = Observable(data["w"][:, iy, :, 1])

q_vals  = data["q"][:, iy, :, :]
n_vals = data["n"][:, iy, :, :]
qv_vals  = data["qv"][:, iy, :, :]
w_vals  = data["w"][:, iy, :, :]

n_eff = max.(n_vals, 0.0)
q_eff = max.(q_vals, 0.0)

p23 = 2.0 / 3.0

tau = similar(n_eff, Float64)

mask = (n_eff .> 0.0) .& (q_eff .> 0.0)

tau .= NaN

tau[mask] .= tref .* mref^p23 .* tau_star .* n_eff[mask].^p23 .* q_eff[mask].^(-p23)

tau_finite = tau[isfinite.(tau)]

tau_min, tau_max = extrema(tau_finite)
println("Computed tau min/max = ", tau_min, " ", tau_max)
println("Median tau = ", median(tau_finite))

tau_obs = Observable(tau[:, iy, :, 1])

#println("q =", q_vals)

q_min  = minimum(skipmissing(filter(isfinite, q_vals)))
q_max  = maximum(skipmissing(filter(isfinite, q_vals)))

qv_min  = minimum(skipmissing(filter(isfinite, qv_vals)))
qv_max  = maximum(skipmissing(filter(isfinite, qv_vals)))

println("q min/max  = ", q_min, " ", q_max)
println("qv min/max  = ", qv_min, " ", qv_max)

n_min = minimum(skipmissing(filter(isfinite, n_vals)))
n_max = maximum(skipmissing(filter(isfinite, n_vals)))

println("n min/max  = ", n_min, " ", n_max)

w_min = minimum(skipmissing(filter(isfinite, w_vals)))
w_max = maximum(skipmissing(filter(isfinite, w_vals)))

fig = Figure(size = (1200, 800))

ax1 = Axis(fig[1,1], xlabel="x (km)", ylabel="z (km)", title="q")
ax2 = Axis(fig[2,1], xlabel="x (km)", ylabel="z (km)", title="n")
ax3 = Axis(fig[3,1], xlabel="x (km)", ylabel="z (km)", title="tau")
ax4 = Axis(fig[4,1], xlabel="x (km)", ylabel="z (km)", title="qv")
ax5 = Axis(fig[5,1], xlabel="x (km)", ylabel="z (km)", title="w")

hm1 = heatmap!(ax1, x, z, q_obs; colormap=:blues, colorrange=(q_min, q_max))
hm2 = heatmap!(ax2, x, z, qv_obs; colormap=:blues, colorrange=(qv_min, qv_max))
hm3 = heatmap!(ax3, x, z, n_obs; colormap=:blues, colorrange=(n_min, n_max))
hm4 = heatmap!(ax4, x, z, tau_obs; colormap=:blues, colorrange=(tau_min, tau_max), nan_color = RGBAf(0,0,0,0))
hm5 = heatmap!(ax5, x, z, w_obs; colormap=:RdBu, colorrange=(w_min, w_max))

Colorbar(fig[1,2], hm1, label="q [kg/kg]")
Colorbar(fig[2,2], hm2, label="qv [kg/kg]")
Colorbar(fig[3,2], hm3, label="n [kg/kg]")
Colorbar(fig[4,2], hm4, label="tau [s]")
Colorbar(fig[5,2], hm5, label="w [m/s]")

tau_frame = similar(q_obs[], Float64)

record(fig,
       output_file,
       1:tmax;
       framerate = 10) do t

    q = data["q"][:, iy, :, t]
    n = data["n"][:, iy, :, t]
    qv = data["qv"][:, iy, :, t]

    q_obs[] = q
    n_obs[] = n
    qv_obs[] = qv
    w_obs[] = data["w"][:, iy, :, t]

    # Effektive Werte
    q_eff = max.(q, 0.0)
    n_eff = max.(n, 0.0)

    # Maske
    mask = (q_eff .> 0.0) .& (n_eff .> 0.0)

    # tau neu berechnen
    tau_frame .= NaN
    tau_frame[mask] .= tref .* mref .* tau_star .*
                        n_eff[mask].^(2/3) .*
                        q_eff[mask].^(-2/3)

    tau_obs[] = tau_frame

    ax1.title = "q at t = $t"
    ax2.title = "qv at t = $t"
    ax3.title = "n at t = $t"
    ax4.title = "tau at t = $t"
    ax5.title = "w at t = $t"
end
close(data)

println("GIF saved to $output_file")