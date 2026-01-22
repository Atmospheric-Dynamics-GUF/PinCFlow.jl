# RUN in PF/pinc/
# julia --project=. examples/visualization/fast_plot.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

# Base font-size (points). Change this to scale all plot text.
FS = 12

data = h5open("/home/dolaptch/PF/pinc/mountain_wave.h5")
#data = h5open("/home/dolaptch/PF/pinc/exp/pincflow_output.h5")
#data = h5open("/home/dolaptch/PF/runs/tjl05/pincflow_output.h5")
#data = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl05/pincflow_output.h5")

#data2s = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5")
data2s = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5")

# Set grid for data1
x = data["x"][:] .* 0.001 
y = data["y"][:] .* 0.001
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Set grid for data2
 x2 = data2s["x2"][:] .* 0.001 
 y2 = data2s["y2"][:] .* 0.001 
 z2 = data2s["z2"][:, :, :] .* 0.001
 x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
 y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

iy = Int(length(data["w"][1, :, 1, 1]) / 2)
iys = 1 
println("iy = ", iy)
tidx = length(data["w"][1, 1, 1, :])
tidx2s = length(data2s["w"][1, 1, 1, :])

fld = data["iaux1"][:, :, :, tidx]
fld2s = data2s["sn"][:, :, :, tidx2s]

println("size fld  ", size(fld))
println("size fld2s", size(fld2s))
println("fld  ", maximum(fld), " ", minimum(fld))
#println("fld2s", maximum(fld2s), " ", minimum(fld2s))

# Plot comparison: data vs datas
figure(figsize=(12, 10))
subplot(211)
contour = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
cbar = colorbar(contour)
cbar.ax.tick_params(labelsize=FS)
title("data: Res", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
ax = gca(); ax.tick_params(labelsize=FS)

subplot(212)
contour2s = pcolormesh(x2[:, iys, :], z2[:, iys, :], fld2s[:,iys,:]; cmap="Blues")
cbar2 = colorbar(contour2s)
cbar2.ax.tick_params(labelsize=FS)
title("data2s: Par (saved)", fontsize=1.5*FS)
xlabel("x (km)", fontsize=1.2*FS); ylabel("z (km)", fontsize=1.2*FS)
ax = gca(); ax.tick_params(labelsize=FS)

tight_layout()
savefig("/home/dolaptch/output/mountain_wave.png")

##############################
# Compute and plot PDFs (probability density estimates) of the field data
##############################

var_name = "\$S\$"
label_name = "saturation ratio $var_name"

#var_name = "\$\\theta'\$ [K]"
#label_name = "potential temperature $var_name"
#fld = data["thetap"][:, :, :, tidx]

#var_name = "\$w\$ [m s\$^{-1}\$]"
#label_name = "vertical velocity $var_name"
#fld = data["w"][:, :, :, tidx]

v1_full = vec(fld)
v2 = vec(fld2s)

# Apply mask from data["n"] at the same time index: remove points where n == 0
nfield = vec(data["n"][:, :, :, tidx])
total_n_v1 = length(v1_full)
keep_mask = nfield .< 100.0

#plot all/ mask only
#v1 = v1_full[keep_mask]
v1 = v1_full
after_mask_n_v1 = length(v1)

# filter out non-finite values
v1 = v1[isfinite.(v1)]
v2 = v2[isfinite.(v2)]
after_finite_n_v1 = length(v1)

# filter out values equal to 0 or 1.5 from fld (v1)
orig_n_v1 = after_finite_n_v1
v1 = v1[(v1 .!= 0.0) .& (v1 .!= 1.5)]
filtered_n_v1 = length(v1)

println("samples fld (total grid points): ", total_n_v1)
println("samples fld (after n<100 mask): ", after_mask_n_v1)
println("samples fld (after non-finite filter): ", after_finite_n_v1)
println("samples fld (after removing 0 and 1.5): ", filtered_n_v1)
println("samples fld2s: ", length(v2), ", mean/std: ", mean(v2), ", ", std(v2))

nbins = 100
figure(figsize=(8, 6))
PyPlot.hist(v1; bins=nbins, density=true, alpha=0.6)

# Fit Gaussian using sample mean and std and overplot the PDF
if length(v1) > 5
	mu = mean(v1)
	sigma = std(v1)
	println("Gaussian fit: mu=", mu, ", sigma=", sigma)
	if sigma > 0
		xg = range(minimum(v1), stop=maximum(v1), length=400)
		yg = 1.0/(sigma*sqrt(2*pi)) .* exp.(-0.5 .* ((xg .- mu)./sigma).^2)
		PyPlot.plot(collect(xg), collect(yg); color="k", linewidth=2, label="Gaussian fit")
	end
end

tight_layout()
legend(fontsize=FS)
title("PDF of \$S(T)\$", fontsize=1.5*FS)
#title("PDF of $var_name in the cloud-free region", fontsize=1.5*FS)
xlabel(label_name, fontsize=1.2*FS)
ylabel("probability density", fontsize=1.2*FS)
ax = gca(); ax.tick_params(labelsize=FS)
tight_layout()
tight_layout()
savefig("/home/dolaptch/output/mountain_wave_pdf.png")

