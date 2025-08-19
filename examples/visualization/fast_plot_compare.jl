# RUN in PF/pinc/
# julia --project=. examples/visualization/fast_plot.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics
using NCDatasets

data = h5open("/home/dolaptch/PF/pinc/test/pincflow_output.h5")
data2 = Dataset("/home/dolaptch/PF/runs/tjl01/pincflow_data_out.nc", "r")

# Set grid.
x = data["x"][:] .* 0.001 .- 10
y = data["y"][:] .* 0.001 .- 10
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]



iy = 1
tidx = size(data["w"][:, :, :, :], 4)  # for a 4D array, 4th dimension is often "time"
#tidx = 1

fld = data["w"][:, :, :, tidx]
#fld = data["n"][:, :, :, tidx]
#fld = data["iaux3"][:, :, :, tidx]

fld2 = data2.group["atmvar"]["w"][:, :, :, tidx]
#fld2 = data2.group["icevar"]["n"][:, :, :, tidx]
#fld2 = data2.group["optvar"]["s"][:, :, :, tidx] 
#fld2 = data2.group["optvar"]["w"][:, :, :, tidx] 
#fld2 = data2.group["optvar"]["t"][:, :, :, tidx] 

println(maximum(fld), " ", minimum(fld))
println(maximum(fld2), " ", minimum(fld2))

subplot(131)
subplots_adjust(wspace=1.6) 
#fld = data["iaux1"][:, :, :, :]
#fld = data["iaux1"][:, :, :, 1]

#contour = contourf( fld[:,1,:]')
#contour = contourf(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
colorbar(contour)

#plot(fld[Int(length(fld[:,1,1])//2),1,:])
#plot(mean(fld[:,1,:], dims=1), label="Mean over x")

subplot(132)

contour2 = pcolormesh(x[:, iy, :], z[:, iy, :], fld2[:,iy,:]; cmap="Blues")
colorbar(contour2)

subplot(133)
contour2 = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]-fld2[:,iy,:]; cmap="Blues")
colorbar(contour2)

savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave.png")

clf()
plot(fld[1,iy,:], label="fld")
plot(fld2[1,iy,:], label="fld2")
legend()
savefig("/home/dolaptch/PF/pinc/examples/visualization/test.png")
