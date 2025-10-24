# RUN in PF/pinc/
# julia --project=. examples/visualization/fast_plot_compare.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics
using NCDatasets

#data = h5open("test/pseudo_incompressible_isothermal_mountain_wave.h5")# canonical test
data = h5open("/home/dolaptch/PF/pinc/test/pincflow_output.h5")
#data = h5open("/home/dolaptch/PF/runs/pjl01/pincflow_output.h5") # v1 in julia
#*data2 = h5open("/home/dolaptch/PF/runs/pjl03/pincflow_output.h5") # v2 in julia change clf according to fortran namelist
#data2 = Dataset("/home/dolaptch/PF/runs/tjl01/pincflow_data_out.nc", "r") # fortran canonical test
data2 = Dataset("/home/dolaptch/PF/runs/tjl04/pincflow_data_out.nc", "r") # 

# Set grid.
x = data["x"][:] .* 0.001 
y = data["y"][:] .* 0.001 
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

x2 = data["x2"][:] .* 0.001 
y2 = data["y2"][:] .* 0.001 
z2 = data["z2"][:, :, :] .* 0.001

# x2 = data["x"][:] .* 0.001
# y2 = data["y"][:] .* 0.001
# z2 = data["z"][:, :, :] .* 0.001

x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

iy = 1
#tidx = size(data["w"][:, :, :, :], 4)  # for a 4D array, 4th dimension is often "time"
tidx = length(data["w"][1, 1, 1, :])  # for a 4D array, 4th dimension is often "time"

#fld = data["rhobar"][:, :, :]
#fld = data["thetabar"][:, :, :]
#fld = data["u"][:, :, :, tidx]
fld = data["sn"][:, :, :, tidx]
#fld = data["wwp"][:, :, :, tidx] #+ data["w"][:, :, :, tidx]
#fld = data["thetap"][:, :, :, tidx]
#fld = data["rhop"][:, :, :, tidx]
#fld = data["n"][:, :, :, tidx]
#fld = data["iaux3"][:, :, :, tidx]
#fld = data["n2"][:, :, :]
#fld0 = data["pip"][:, :, :, 1]
#fld = data["pip"][:, :, :, tidx]
#fld = data["thetabar"][:, :, :]
println("size fld ", size(fld))

#fld2 = data2.group["atmvar"]["rhobar"][:]
#fld2 = data2.group["atmvar"]["thetabar"][:]
#fld2 = data2.group["atmvar"]["rhop"][:, :, :, tidx]
#fld2 = data2.group["atmvar"]["thetap"][:, :, :, tidx]
#fld2 = data2.group["atmvar"]["u"][:, :, :, tidx]
fld2 = data2.group["icevar"]["sn"][:, :, :, tidx]
#fld2 = data2.group["icevar"]["sqv"][:, :, :, tidx]
#fld2 = data2.group["optvar"]["w"][:, :, :, tidx]
#fld20 = data2.group["atmvar"]["pip"][:, :, :, 1]
#fld2 = data2.group["atmvar"]["pip"][:, :, :, tidx]
#fld2 = data2.group["atmvar"]["thetabar"][:]
#fld2 = data2.group["atmvar"]["N2"][:, :, :]
#fld2 = data2.group["optvar"]["s"][:, :, :, tidx] 
#fld2 = data2.group["optvar"]["w"][:, :, :, tidx] 
#fld2 = data2.group["optvar"]["t"][:, :, :, tidx] 

println("size fld2 ", size(fld2))
println(".jl  ", maximum(fld), " ", minimum(fld))
println(".f90 ", maximum(fld2), " ", minimum(fld2))

subplot(131)
subplots_adjust(wspace=1.6) 
#fld = data["iaux1"][:, :, :, :]
#fld = data["iaux1"][:, :, :, 1]

yaxis_limits = (7, 9.5)
ylim(yaxis_limits)
#contour = contourf( fld[:,1,:]')
#contour = contourf(x2[:, iy, :], z2[:, iy, :], fld[:,iy,:]; cmap="RdBu")
contour = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld[:,iy,:]; cmap="Blues")
colorbar(contour)
title(".jl")
#plot(fld[Int(length(fld[:,1,1])//2),1,:])

#plot(mean(fld[:,1,:], dims=1), label="Mean over x")

subplot(132)

#only if thetabar/rhobar
#fld2 = repeat(fld2, 1, size(fld)[1])'
#fld2 = reshape(fld2, size(fld2)[1], 1, size(fld2)[2])

#contour2 = contourf( fld2[:,1,:]')
contour2 = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
#contour2 = contourf(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="RdBu")
colorbar(contour2)
ylim(yaxis_limits)
title(".f90")

subplot(133)
#contour2 = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]-fld2[:,iy,:]; cmap="Blues")
contour2 = contourf(x2[:, iy, :], z2[:, iy, :], fld[:,iy,:]-fld2[:,iy,:]; cmap="RdBu")
colorbar(contour2)
ylim(yaxis_limits)
title("Difference")
savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave.png")


clf()
ix = 1
iz=15
plot(fld[1,iy,:], label="fld")

#plot(fld2[ix,iy,:], label="fld2")

#plot(fld[:,iy,iz], label="fld")
#plot(fld2[:,iy,iz], label="fld2")
#plot(fld[:,1,10]-fld2[:,1,10], label="fld2")

clf()
ix = 1
iz = 160
ix = 1
#xvals = 1:length(fld[:,iy,iz])
xvals = 1:length(fld[ix,iy,:])

fig, ax1 = subplots()
ax1.plot(xvals, fld[ix,iy,:], "b--", label="fld")
ax1.plot(xvals, fld2[ix,iy,:], "k--", label="fld")
#ax1.set_ylabel("fld", color="b")
ax1.set_ylabel("fld", color="b")

ax2 = ax1.twinx()
#ax2.plot(xvals, abs.(fld[ix,iy,:]) - abs.(fld2[ix,iy,:]), "r-", label="diff")
#ax2.set_ylabel("diff", color="r")

#legend()
#savefig("/home/dolaptch/PF/pinc/examples/visualization/test_twin_axis.png")
#plot(fld[:,1,10], ls = "--", label="fld")

legend()
savefig("/home/dolaptch/PF/pinc/examples/visualization/test.png")
