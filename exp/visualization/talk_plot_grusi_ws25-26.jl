# RUN in PF/pinc/
# julia --project=. examples/visualization/talk_plot_grusi_ws25-26.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot
using Statistics

data = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl05/pincflow_output.h5")  #LES
data2 = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl06/pincflow_output.h5") # RT + SGS Ice

# Set grid.
x = data["x"][:] .* 0.001
y = data["y"][:] .* 0.001
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

x2 = data2["x2"][:] .* 0.001
y2 = data2["y2"][:] .* 0.001
z2 = data2["z2"][:, :, :] .* 0.001

x2 = [xi for xi in x2, iy in 1:size(z2)[2], iz in 1:size(z2)[3]]
y2 = [yi for ix in 1:size(z2)[1], yi in y2, iz in 1:size(z2)[3]]

iy = 1
tidx = length(data["w"][1, 1, 1, :])  # for a 4D array, 4th dimension is often "time"

fld = data["n"][:, :, :, tidx]
fld2 = data2["sn"][:, :, :, tidx]
#fld2 = data2["wwp"][:, :, :, tidx] #+ data["w"][:, :, :, tidx]

println("size fld ", size(fld))
println("size fld2 ", size(fld2))
println("fld  ", maximum(fld), " ", minimum(fld))
println("fld2 ", maximum(fld2), " ", minimum(fld2))

# ==========================
# Plot ice number
# ==========================

figure()
# Figure-level title spanning both subplots
subplots_adjust(top = 0.90)
suptitle("Ice Number Concentration")

ax1 = subplot(121)

yaxis_limits = (0.0, 12.0)
# Use the same color scale for both panels
vmin = min(minimum(fld), minimum(fld2))
vmax = max(maximum(fld), maximum(fld2))
ylim(yaxis_limits)

im1 = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:, iy, :]; cmap = "Blues", vmin = vmin, vmax = vmax)
contour(x[:, iy, :], z[:, iy, :], fld[:, iy, :], levels = [0.0], colors = "k", linewidths = 0.5)
title("Res");
ylabel("z [km]");
xlabel("x [km]")

ax2 = subplot(122)
h2 = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2[:, iy, :]; cmap = "Blues", vmin = vmin, vmax = vmax)
#h2 = contourf(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
contour(x2[:, iy, :], z2[:, iy, :], fld2[:, iy, :], levels = [0.0], colors = "k", linewidths = 0.5)
#clim(vlow, vup)
ylim(yaxis_limits)
title("SGS Ice Par GW");
xlabel("x [km]");

# Single, shared colorbar for both panels
colorbar(h2, ax = [ax1, ax2])

savefig("/home/dolaptch/PF/pinc/examples/visualization/nice_res_par.png")

# ==========================
# Histogram comparison (LES vs RT) for slices fld[:, iy, :] and fld2[:, iy, :]
# ==========================

figure()

# Figure-level title spanning both subplots
subplots_adjust(top = 0.90)
suptitle("Ice Number Histogram")

# Flatten values along the selected j-index and keep finite values only
vals1 = vec(fld[:, iy, :])
vals2 = vec(fld2[:, iy, :])
finite1 = vals1[isfinite.(vals1)]
finite2 = vals2[isfinite.(vals2)]

# Shared linear bins for comparable histograms
if !isempty(finite1) && !isempty(finite2)
	low1 = 1.0e2 # critical value to define a cloud
	low2 = 1.0e2 # critical value to define a cloud
	low = min(minimum(finite1), minimum(finite2))
	high = max(maximum(finite1), maximum(finite2))
	high1 = maximum(finite1)
	high2 = maximum(finite2)
	if low == high
		# Expand a degenerate range slightly to avoid errors
		δ = max(abs(high), 1.0) * 1e-6
		low -= δ;
		high += δ
	end
	edges = collect(range(low, high; length = 11))
	edges1 = collect(range(low1, high1; length = 11))
	edges2 = collect(range(low2, high2; length = 11))
	println("Histogram edges2: ", edges2)
	subplot(1, 2, 1)
	hist(finite1; bins = edges1, alpha = 0.5, density = true, label = "LES", color = "C0")
	xlabel("n [1/kg]");
	title("Res.")
	subplot(1, 2, 2)
	hist(finite2; bins = edges2, alpha = 0.5, density = true, label = "RT", color = "C1")
	xlabel("n [1/kg]");
	title("SGS Ice Par GW.")
	tight_layout()
	savefig("/home/dolaptch/PF/pinc/examples/visualization/nice_res_par_hist.png")
else
	println("No finite values found for histogram plots.")
end

# ==========================
# Function: coarse-grid fraction of exceedances
# ==========================

function coarse_fraction(X::AbstractArray, Z::AbstractArray, F::AbstractArray; thresh::Real = 0.0, nxc::Int = 40, nzc::Int = 30)
	xmin, xmax = minimum(X), maximum(X)
	zmin, zmax = minimum(Z), maximum(Z)
	if xmin == xmax
		δ = max(abs(xmax), 1.0) * 1e-6;
		xmin -= δ;
		xmax += δ
	end
	if zmin == zmax
		δ = max(abs(zmax), 1.0) * 1e-6;
		zmin -= δ;
		zmax += δ
	end
	xedges = collect(range(xmin, xmax; length = nxc+1))
	zedges = collect(range(zmin, zmax; length = nzc+1))

	frac = zeros(Float64, nxc, nzc)
	for ii in 1:nxc
		xlo, xhi = xedges[ii], xedges[ii+1]
		xmask = (X .>= xlo) .& ((ii < nxc) ? (X .< xhi) : (X .<= xhi))
		for jj in 1:nzc
			zlo, zhi = zedges[jj], zedges[jj+1]
			zmask = (Z .>= zlo) .& ((jj < nzc) ? (Z .< zhi) : (Z .<= zhi))
			mask = xmask .& zmask
			total = count(mask)
			if total == 0
				frac[ii, jj] = NaN
			else
				exceed = count((F .> thresh) .& mask)
				frac[ii, jj] = exceed / total
			end
		end
	end

	xc = (xedges[1:(end-1)] .+ xedges[2:end]) ./ 2
	zc = (zedges[1:(end-1)] .+ zedges[2:end]) ./ 2
	return xc, zc, frac
end

function fraction_pdf(frac::AbstractArray; nb::Int = 50)
	vals = vec(frac)
	v = vals[isfinite.(vals) .& (vals .>= 0) .& (vals .<= 1)]
	edges = collect(range(0.0, 1.0; length = nb+1))
	counts = zeros(Int, nb)
	for val in v
		bi = searchsortedlast(edges, val) - 1
		bi = clamp(bi, 1, nb)
		counts[bi] += 1
	end
	binw = 1.0 / nb
	total = sum(counts)
	pdf = total == 0 ? zeros(Float64, nb) : counts ./ (total * binw)
	centers = (edges[1:(end-1)] .+ edges[2:end]) ./ 2
	return centers, pdf
end

# Parameters (adjust as needed or wire to ARGS)
thresh = 0.0                  # threshold for exceedance
nxc, nzc = 8, 8            # coarse grid resolution in x and z

# Compute coarse fractions for both LES and RT using the function and compare PDFs
xc_rt, zc_rt, frac_rt = coarse_fraction(x2[:, iy, :], z2[:, iy, :], fld2[:, iy, :]; thresh = thresh, nxc = nxc, nzc = nzc)
xc_ls, zc_ls, frac_ls = coarse_fraction(x[:, iy, :], z[:, iy, :], fld[:, iy, :]; thresh = thresh, nxc = nxc, nzc = nzc)

figure()
subplots_adjust(top = 0.90)
suptitle("Cloud Fraction (5 x 2 km)")

# Coarse cloud fraction maps (Res vs Par) with shared colorbar
ax1 = subplot(121)
pc1 = pcolormesh(xc_ls, zc_ls, frac_ls'; cmap = "Blues", vmin = 0.0, vmax = 1.0)
title("Res")
xlabel("x (km)");
ylabel("z (km)")

ax2 = subplot(122)
pc2 = pcolormesh(xc_rt, zc_rt, frac_rt'; cmap = "Blues", vmin = 0.0, vmax = 1.0)
colorbar(pc2, ax = [ax1, ax2])
title("SGS Ice Par GW")
xlabel("x (km)")
savefig("/home/dolaptch/PF/pinc/examples/visualization/ccf_res_par.png")

c_rt, pdf_rt = fraction_pdf(frac_rt; nb = 50)
c_ls, pdf_ls = fraction_pdf(frac_ls; nb = 50)

mask_rt = pdf_rt .>= 0
mask_ls = pdf_ls .>= 0

figure()
plot(c_ls[mask_ls], pdf_ls[mask_ls]; label = "Res", color = "C0")
plot(c_rt[mask_rt], pdf_rt[mask_rt]; label = "Par", color = "C1")
xlabel("cloud-cover fraction")
#ylabel("pdf")
title("PDF of cloud-cover fractions")
legend()
grid(true)
tight_layout()
savefig("/home/dolaptch/PF/pinc/examples/visualization/ccf_res_par_pdf.png")

# ==========================
# compare with fully parameterized ice-microphysics
# ==========================

data3 = h5open("/home/dolaptch/PF/runs_save/GWI_grusi_WS_25-26/tjl08/pincflow_output.h5") # RT + cloud cover

tidx = 3
println("Time index for noSGS changed to ", tidx)
sleep(3)
clc = data3["clc"][:, :, :, tidx]

figure()
subplots_adjust(top = 0.90)
suptitle("Cloud Fraction (5 x 2 km)")

# Coarse cloud fraction maps (Res vs Par) with shared colorbar
ax1 = subplot(131)
pc1 = pcolormesh(xc_ls, zc_ls, frac_ls'; cmap = "Blues", vmin = 0.0, vmax = 1.0)
title("Res")
xlabel("x (km)");
ylabel("z (km)")

ax2 = subplot(132)
pc2 = pcolormesh(xc_rt, zc_rt, frac_rt'; cmap = "Blues", vmin = 0.0, vmax = 1.0)
title("SGS Ice Par GW")
xlabel("x (km)")
gca().set_yticklabels([])

ax3 = subplot(133)
pc3 = pcolormesh(xc_rt, zc_rt, clc[:, iy, :]'; cmap="Blues", vmin=0.0, vmax=1.0)
colorbar(pc3, ax = [ax1, ax2, ax3])
title("Par Ice+GW")
xlabel("x (km)")
gca().set_yticklabels([])
savefig("/home/dolaptch/PF/pinc/examples/visualization/ccf_res_par_nosgs.png")
