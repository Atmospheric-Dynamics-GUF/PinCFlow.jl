# RUN in PF/pinc/
# julia --project=. examples/visualization/fast_plot_compare.jl
# view with: eog examples/visualization/mountain_wave.png &
#
using HDF5
#using LaTeXStrings
using PyPlot 
using Statistics
using NCDatasets

data = h5open("/home/dolaptch/PF/runs/tjl05/pincflow_output.h5")  #LES
data2 = h5open("/home/dolaptch/PF/runs/tjl06/pincflow_output.h5") # RT
#data2 = h5open("/home/dolaptch/PF/pinc/test/pincflow_output.h5") # RT

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
tidx = 30 #length(data["w"][1, 1, 1, :])  # for a 4D array, 4th dimension is often "time"

#fld = data["rhobar"][:, :, :]
#fld = data["thetabar"][:, :, :]
fld = data["n"][:, :, :, tidx]
fld2 = data2["sn"][:, :, :, tidx]
#fld2 = data2["wwp"][:, :, :, tidx] #+ data["w"][:, :, :, tidx]
#fld = data["thetap"][:, :, :, tidx]
#fld = data["rhop"][:, :, :, tidx]
#fld = data["n"][:, :, :, tidx]
#fld2 = data2["iaux1"][:, :, :, tidx]
#fld = data["n2"][:, :, :]
#fld0 = data["pip"][:, :, :, 1]
#fld = data["pip"][:, :, :, tidx]
#fld = data["thetabar"][:, :, :]
println("size fld ", size(fld))
println("size fld2 ", size(fld2))
println("fld  ", maximum(fld), " ", minimum(fld))
println("fld2 ", maximum(fld2), " ", minimum(fld2))

subplot(131)
subplots_adjust(wspace=1.6) 

#yaxis_limits = (7, 9.5)
#yaxis_limits = (7.5, 10.5)
yaxis_limits = (.0, 12.0)
vlow = minimum(fld); vup = maximum(fld)
#vlow = -0.7
#vup = 0.7
ylim(yaxis_limits)

#contour = contourf(x2[:, iy, :], z2[:, iy, :], fld[:,iy,:]; cmap="RdBu")
h1 = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]; cmap="Blues")
#clim(vlow, vup)  
colorbar(h1)
contour(x[:, iy, :], z[:, iy, :], fld[:,iy,:], levels=[0.0], colors="k", linewidths=0.5)
title("LES")

subplot(132)
h2 = pcolormesh(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
#h2 = contourf(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:]; cmap="Blues")
colorbar(h2)
contour(x2[:, iy, :], z2[:, iy, :], fld2[:,iy,:], levels=[0.0], colors="k", linewidths=0.5)
#clim(vlow, vup)
ylim(yaxis_limits)
title("RT")

subplot(133)
#contour2 = pcolormesh(x[:, iy, :], z[:, iy, :], fld[:,iy,:]-fld2[:,iy,:]; cmap="Blues")
#contour2 = contourf(x2[:, iy, :], z2[:, iy, :], fld[:,iy,:]-fld2[:,iy,:]; cmap="RdBu")
#colorbar(contour2)
#ylim(yaxis_limits)
title("Difference")
savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave.png")


# clf()
# iz = 120
# plot(x[:, iy, iz], fld[:, iy, iz], "-", label="LES")
# plot(x2[:, iy, round(Int, (iz-1)*4+2)], fld2[:, iy, round(Int, (iz-1)*4+2)], "r-", label="RT")
# legend()
# savefig("/home/dolaptch/PF/pinc/examples/visualization/test.png")

# ==========================
# Histogram comparison (LES vs RT) for slices fld[:, iy, :] and fld2[:, iy, :]
# ==========================

figure()

# Flatten values along the selected j-index and keep finite values only
vals1 = vec(fld[:, iy, :])
vals2 = vec(fld2[:, iy, :])
finite1 = vals1[isfinite.(vals1)]
finite2 = vals2[isfinite.(vals2)]

# Shared linear bins for comparable histograms
if !isempty(finite1) && !isempty(finite2)
    low1 = 1.0e2 # critical value to define a cloud
    low2 = 1.0e6 # critical value to define a cloud
	low = min(minimum(finite1), minimum(finite2))
    high = max(maximum(finite1), maximum(finite2))
	high1 = maximum(finite1)
    high2 = maximum(finite2)
	if low == high
		# Expand a degenerate range slightly to avoid errors
		δ = max(abs(high), 1.0) * 1e-6
		low -= δ; high += δ
	end
	edges = collect(range(low, high; length=11))
    edges1 = collect(range(low1, high1; length=11))
    edges2 = collect(range(low2, high2; length=11))
    println("Histogram edges2: ", edges2) 
	subplot(1, 2, 1)
	hist(finite1; bins=edges1, alpha=0.5, density=true, label="LES", color="C0")
   subplot(1, 2, 2) 
	hist(finite2; bins=edges2, alpha=0.5, density=true, label="RT", color="C1")
	xlabel("Value"); ylabel("PDF"); title("Histogram (linear)")
	legend()

	# # Log-x histogram (only positive values)
	# subplot(1, 2, 2)
	# pos1 = finite1[finite1 .> 0]
	# pos2 = finite2[finite2 .> 0]
	# if !isempty(pos1) && !isempty(pos2)
	# 	lowp = min(minimum(pos1), minimum(pos2))
	# 	highp = max(maximum(pos1), maximum(pos2))
	# 	if lowp < highp
	# 		edges_log = 10 .^ range(log10(lowp), log10(highp); length=51)
	# 		hist(pos1; bins=edges_log, alpha=0.5, density=true, label="LES", color="C0")
	# 		hist(pos2; bins=edges_log, alpha=0.5, density=true, label="RT", color="C1")
	# 		xscale("log")
	# 		xlabel("Value (log)"); ylabel("PDF"); title("Histogram (log-x)")
	# 		legend()
	# 	else
	# 		title("Histogram (log-x): insufficient range")
	# 	end
	# else
	# 	title("Histogram (log-x): no positive values")
	# end

	tight_layout()
	savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave_hist.png")
else
	println("No finite values found for histogram plots.")
end

	# ==========================
	# Coarse-grid fraction of exceedances for RT slice fld2[:, iy, :]
	# ==========================

	

	# # Extract 2D fields for the selected slice
	# X = x2[:, iy, :]
	# Z = z2[:, iy, :]
	# F = fld2[:, iy, :]

	# # Compute bin edges from coordinate ranges
	# xmin, xmax = minimum(X), maximum(X)
	# zmin, zmax = minimum(Z), maximum(Z)
	# if xmin == xmax
	# 	δ = max(abs(xmax), 1.0) * 1e-6; xmin -= δ; xmax += δ
	# end
	# if zmin == zmax
	# 	δ = max(abs(zmax), 1.0) * 1e-6; zmin -= δ; zmax += δ
	# end
	# xedges = collect(range(xmin, xmax; length=nxc+1))
	# zedges = collect(range(zmin, zmax; length=nzc+1))

	# # Initialize fraction grid
	# frac = zeros(Float64, nxc, nzc)

	# # Precompute masks per coarse cell and fill fractions
	# for ii in 1:nxc
	# 	xlo, xhi = xedges[ii], xedges[ii+1]
	# 	# include right edge only for last bin to cover the max
	# 	xmask = (X .>= xlo) .& ((ii < nxc) ? (X .< xhi) : (X .<= xhi))
	# 	for jj in 1:nzc
	# 		zlo, zhi = zedges[jj], zedges[jj+1]
	# 		zmask = (Z .>= zlo) .& ((jj < nzc) ? (Z .< zhi) : (Z .<= zhi))
	# 		mask = xmask .& zmask
	# 		total = count(mask)
	# 		if total == 0
	# 			frac[ii, jj] = NaN
	# 		else
	# 			exceed = count( (F .> thresh) .& mask )
	# 			frac[ii, jj] = exceed / total
	# 		end
	# 	end
	# end

	# # Plot the coarse fraction map using bin centers
	# xc = (xedges[1:end-1] .+ xedges[2:end]) ./ 2
	# zc = (zedges[1:end-1] .+ zedges[2:end]) ./ 2

	# figure()
	# pc = pcolormesh(xc, zc, frac'; cmap="Blues", vmin=0.0, vmax=1.0)
	# colorbar(pc)
	# xlabel("x (km)"); ylabel("z (km)"); title("Fraction > $(thresh) in coarse cells (RT)")
	# savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave_fraction.png")

	# # Histogram of fraction values and log-PDF plot
	# fvals = vec(frac)
	# fvalid = fvals[isfinite.(fvals) .& (fvals .>= 0.0) .& (fvals .<= 1.0)]
	# if !isempty(fvalid)
	# 	nb = 50
	# 	edgesf = collect(range(0.0, 1.0; length=nb+1))
	# 	counts = zeros(Int, nb)
	# 	for v in fvalid
	# 		bi = searchsortedlast(edgesf, v) - 1
	# 		if 1 <= bi <= nb
	# 			counts[bi] += 1
	# 		end
	# 	end
	# 	binw = 1.0 / nb
	# 	pdf = counts ./ (sum(counts) * binw)
	# 	binc = (edgesf[1:end-1] .+ edgesf[2:end]) ./ 2
	# 	#maskp = pdf .> 0
    #     maskp = (pdf .>= 0)
	# 	figure()
	# 	plot(binc[maskp], pdf[maskp]; color="C3")
    #     #plot(binc[maskp], log.(pdf[maskp]); color="C3")
	# 	xlabel("fraction (> $(thresh))")
	# 	ylabel("log(pdf)")
	# 	title("Log PDF of coarse-cell fraction (RT)")
	# 	grid(true)
	# 	tight_layout()
	# 	savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave_fraction_logpdf.png")
	# else
	# 	println("No valid fraction values for histogram/log-PDF.")
	# end

		# ==========================
		# Function: coarse-grid fraction of exceedances
		# ==========================

		function coarse_fraction(X::AbstractArray, Z::AbstractArray, F::AbstractArray; thresh::Real=0.0, nxc::Int=40, nzc::Int=30)
			xmin, xmax = minimum(X), maximum(X)
			zmin, zmax = minimum(Z), maximum(Z)
			if xmin == xmax
				δ = max(abs(xmax), 1.0) * 1e-6; xmin -= δ; xmax += δ
			end
			if zmin == zmax
				δ = max(abs(zmax), 1.0) * 1e-6; zmin -= δ; zmax += δ
			end
			xedges = collect(range(xmin, xmax; length=nxc+1))
			zedges = collect(range(zmin, zmax; length=nzc+1))

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

			xc = (xedges[1:end-1] .+ xedges[2:end]) ./ 2
			zc = (zedges[1:end-1] .+ zedges[2:end]) ./ 2
			return xc, zc, frac
		end

		function fraction_pdf(frac::AbstractArray; nb::Int=50)
			vals = vec(frac)
			v = vals[isfinite.(vals) .& (vals .>= 0) .& (vals .<= 1)]
			edges = collect(range(0.0, 1.0; length=nb+1))
			counts = zeros(Int, nb)
			for val in v
				bi = searchsortedlast(edges, val) - 1
				bi = clamp(bi, 1, nb)
				counts[bi] += 1
			end
			binw = 1.0 / nb
			total = sum(counts)
			pdf = total == 0 ? zeros(Float64, nb) : counts ./ (total * binw)
			centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
			return centers, pdf
		end

        # Parameters (adjust as needed or wire to ARGS)
	    thresh = 0.0                  # threshold for exceedance
	    nxc, nzc = 8, 8            # coarse grid resolution in x and z

		# Compute coarse fractions for both LES and RT using the function and compare PDFs
		xc_rt, zc_rt, frac_rt = coarse_fraction(x2[:, iy, :], z2[:, iy, :], fld2[:, iy, :]; thresh=thresh, nxc=nxc, nzc=nzc)
		xc_ls, zc_ls, frac_ls = coarse_fraction(x[:, iy, :],  z[:, iy, :],  fld[:,  iy, :];  thresh=thresh, nxc=nxc, nzc=nzc)

		c_rt, pdf_rt = fraction_pdf(frac_rt; nb=50)
		c_ls, pdf_ls = fraction_pdf(frac_ls; nb=50)

		mask_rt = pdf_rt .>= 0
		mask_ls = pdf_ls .>= 0

		figure()
		plot(c_ls[mask_ls], pdf_ls[mask_ls]; label="LES", color="C0")
		plot(c_rt[mask_rt], pdf_rt[mask_rt]; label="RT", color="C1")
		xlabel("fraction (> $(thresh))")
		ylabel("pdf")
		title("PDF comparison of cloudcover fractions")
		legend()
		grid(true)
		tight_layout()
		savefig("/home/dolaptch/PF/pinc/examples/visualization/mountain_wave_fraction_pdf_compare.png")