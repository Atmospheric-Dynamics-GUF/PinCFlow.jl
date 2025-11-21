# RUN in PF/pinc/
# julia --project=. examples/visualization/fast_plot.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

data = h5open("/home/dolaptch/PF/runs/tjl05/pincflow_output.h5")
datas = h5open("/home/dolaptch/PF/runs/sjl05/pincflow_output.h5") # shift domain to -20,20 km to be consistent with merged run
#datas = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl05/pincflow_output.h5")

data2 = h5open("/home/dolaptch/PF/runs/tjl06/pincflow_output.h5")
data2s = h5open("/home/dolaptch/PF/runs/sjl06/pincflow_output.h5") # shift domain to -20,20 km to be consistent with merged run
#data2s = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5")

# Set grid for data1
x = data["x"][:] .* 0.001 
y = data["y"][:] .* 0.001
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Set grid for data2
x2 = data2["x2"][:] .* 0.001 
y2 = data2["y2"][:] .* 0.001 
z2 = data2["z2"][:, :, :] .* 0.001
x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

iy = 1
tidx = length(data["w"][1, 1, 1, :])
tidxs = tidx
tidx2 = length(data2["w"][1, 1, 1, :])
tidx2s = tidx2 

fld = data["n"][:, :, :, tidx]
flds = datas["n"][:, :, :, tidxs]
fld2 = data2["sn"][:, :, :, tidx2]
fld2s = data2s["sn"][:, :, :, tidx2s]

println("size fld  ", size(fld))
println("size flds ", size(flds))
println("size fld2  ", size(fld2))
println("size fld2s", size(fld2s))
println("fld  ", maximum(fld), " ", minimum(fld))
println("flds ", maximum(flds), " ", minimum(flds))
println("fld2 ", maximum(fld2), " ", minimum(fld2))
println("fld2s", maximum(fld2s), " ", minimum(fld2s))

# Plot comparison: data vs datas
figure(figsize=(12, 10))
subplot(221)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
colorbar(contour)
title("data: Res")
xlabel("x (km)"); ylabel("z (km)")

subplot(222)
contours = pcolormesh(x[:, iy, :], z[:, iy, :], flds[:,iy,:]; cmap="Blues")
colorbar(contours)
title("datas: Res (saved)")
xlabel("x (km)"); ylabel("z (km)")

# Plot comparison: data2 vs data2s
subplot(223)
contour2 = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
colorbar(contour2)
title("data2: Par")
xlabel("x (km)"); ylabel("z (km)")

subplot(224)
contour2s = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2s[:,iy,:]; cmap="Blues")
colorbar(contour2s)
title("data2s: Par (saved)")
xlabel("x (km)"); ylabel("z (km)")

tight_layout()
savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave.png")