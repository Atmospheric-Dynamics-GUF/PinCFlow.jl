using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

data = h5open("/home/dolaptch/PF/pinc/test/pincflow_output.h5")

# Set grid.
x = data["x"][:] .* 0.001 .- 10
y = data["y"][:] .* 0.001 .- 10
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

iy = 1
#fld = data["iaux1"][:, :, :, :]
#println(maximum(fld))
#fld = data["iaux1"][:, :, :, 1]
#fld = data["w"][:, :, :, 1]
fld = data["n"][:, :, :, end]
#contour = contourf( fld[:,1,:]')
contour = contourf(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
colorbar(contour)
#plot(fld[Int(length(fld[:,1,1])//2),1,:])
#plot(mean(fld[:,1,:], dims=1), label="Mean over x")

savefig("/home/dolaptch/PF/pinc/examples/mountain_wave.png")