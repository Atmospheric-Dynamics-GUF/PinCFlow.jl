# RUN in PF/pinc/
# source PF/venv_lev/bin/activate  
# julia --project=exp/ exp/visualization/fast_plot_everywhere.jl
# view with: eog examples/visualization/mountain_wave.png &

using HDF5
#using LaTeXStrings
using PyPlot
using Statistics
using Sockets

# Determine host and user (Julia version)
# runName may be provided as first ARGS entry, otherwise use a default
runName = length(ARGS) >= 1 ? ARGS[1] : "default_run"
host_name = Sockets.gethostname()
user_name = get(ENV, "USER", get(ENV, "LOGNAME", "unknown"))

runName1 = "tjl10"   # Resolved simulation
runName2 = "tjl11"   # RT simulation

file_name = "ice_mountain_wave.h5"

if occursin("levante", host_name)
  # Levante cluster
  dirHome = "/home/b/$user_name/PF/pinc"
  dirScratch = "/scratch/b/$user_name/PF/runs"
  dirSaveCode = "/home/b/$user_name/PF/code_runs"
elseif occursin("login", host_name)
  # Goethe cluster (login node)
  dirHome = "/home/atmodynamics/$user_name/PF/pinc"
  dirScratch = "/scratch/atmodynamics/$user_name/PF/runs"
  dirSaveCode = "/home/atmodynamics/$user_name/PF/code_runs"
else
  # Local machine or unknown
  dirHome = "/home/dolaptch/PF/pinc"
  dirScratch = "/home/dolaptch/PF/runs"
  dirSaveCode = "/home/dolaptch/PF/code_runs"
end

data  = h5open("$dirScratch/$runName1/$file_name")
data2 = h5open("$dirScratch/$runName2/$file_name")

# Set grid for data1
x = data["x"][:] .* 0.001 
y = data["y"][:] .* 0.001
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Set grid for data2
# x2 = data2["x2"][:] .* 0.001 
# y2 = data2["y2"][:] .* 0.001 
# z2 = data2["z2"][:, :, :] .* 0.001
# x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
# y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

x2 = data2["x"][:] .* 0.001 
y2 = data2["y"][:] .* 0.001 
z2 = data2["z"][:, :, :] .* 0.001
x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

iy = 1
tidx = length(data["w"][1, 1, 1, :])
tidxs = tidx
tidx2 = length(data2["w"][1, 1, 1, :])
tidx2s = tidx2 

fld = data["w"][:, :, :, tidx]
fld2 = data2["w"][:, :, :, tidx2]

println("size fld  ", size(fld))
println("size fld2  ", size(fld2))
println("fld  ", maximum(fld), " ", minimum(fld))
println("fld2", maximum(fld2), " ", minimum(fld2))

# Plot comparison: data vs datas
figure(figsize=(12, 10))
subplot(221)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
colorbar(contour)
title("data: Res")
xlabel("x (km)"); ylabel("z (km)")

subplot(222)
contours = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
colorbar(contours)
title("datas: Res (saved)")
xlabel("x (km)"); ylabel("z (km)")

# # Plot comparison: data2 vs data2s
# subplot(223)
# contour2 = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
# colorbar(contour2)
# title("data2: Par")
# xlabel("x (km)"); ylabel("z (km)")

# subplot(224)
# contour2s = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2s[:,iy,:]; cmap="Blues")
# colorbar(contour2s)
# title("data2s: Par (saved)")
# xlabel("x (km)"); ylabel("z (km)")

tight_layout()
savefig("$dirHome/exp/visualization/mountain_wave.png")