include("style.jl")
include("tools.jl")

using LaTeXStrings

# Set paths.
data_path = "../"
reference_path = "../"

println("data_path =", data_path)
println("reference_path =", reference_path)

# Import data.
data = ModelOutput(data_path * "2DWP90_msgwam")
reference = ModelOutput(reference_path * "2DWP90_msgwam")

println("Time: ", data.variables["t"][end], " s")

# Set grid.
x = data.variables["x"][:] .* 0.001
z = data.variables["z"][15:45] .* 0.001
x = [xi for xi in x, j in 1:size(z)[1]]
z = [zj for i in 1:size(x)[1], zj in z]

# Get large-scale tracer distribution change (tls)
# and zonal, meridional, and vertical wind (uls, vls, wls)
tls = data.groups["atmvar"].variables["tmrd"][:, 1, 15:45, end]
uls = data.groups["atmvar"].variables["u"][:, 1, 15:45, end]
vls = data.groups["atmvar"].variables["v"][:, 1, 15:45, end]
wls = data.groups["atmvar"].variables["w"][:, 1, 15:45, end]

# Get reference data.
tlsr = reference.groups["atmvar"].variables["tmrd"][:, 1, 15:45, end]
ulsr = reference.groups["atmvar"].variables["u"][:, 1, 15:45, end]
vlsr = reference.groups["atmvar"].variables["v"][:, 1, 15:45, end]
wlsr = reference.groups["atmvar"].variables["w"][:, 1, 15:45, end]

# Compute differences.
deltat = tls .- tlsr
deltau = uls .- ulsr
deltav = vls .- vlsr
deltaw = wls .- wlsr

# Create figure.
figure(; figsize = (8, 6))

# Plot the zonal wind.
(levels, colormap) = symmetric_contours(minimum(uls), maximum(uls))
subplot(221)
zonal_wind = contourf(x, z, uls; levels = levels, cmap = colormap)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"Large-scale zonal wind $u\,\left[\mathrm{m\,s^{-1}}\right]$")
colorbar(zonal_wind)

# Plot the meridional wind.
(levels, colormap) = symmetric_contours(minimum(vls), maximum(vls))
subplot(222)
meridional_wind = contourf(x, z, vls; levels = levels, cmap = colormap)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"Large-scale meridional wind $v\,\left[\mathrm{m\,s^{-1}}\right]$")
colorbar(meridional_wind)

# Plot the vertical wind.
(levels, colormap) = symmetric_contours(minimum(wls), maximum(wls))
subplot(223)
vertical_wind = contourf(x, z, wls; levels = levels, cmap = colormap)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"Large-scale vertical wind $w\,\left[\mathrm{m\,s^{-1}}\right]$")
colorbar(vertical_wind)

# Plot the tracer.
(levels, colormap) = symmetric_contours(minimum(tls), maximum(tls))
subplot(224)
tracer = contourf(x, z, tls; levels = levels, cmap = colormap)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title("Large-scale tracer")
colorbar(tracer)

# Save the figure.
savefig(data_path * "results/2DWP90_msgwam.png")

# Create difference plots.
if data_path != reference_path
  clf()

  # Plot the zonal-wind-differences.
  (levels, colormap) = symmetric_contours(minimum(deltauls), maximum(deltauls))
  subplot(221)
  zonal_wind = contourf(x, z, deltauls; levels = levels, cmap = colormap)
  xlabel(L"x\,\left[\mathrm{km}\right]")
  ylabel(L"z\,\left[\mathrm{km}\right]")
  title(L"Large-scale zonal wind $u\,\left[\mathrm{m\,s^{-1}}\right]$")
  colorbar(zonal_wind)

  # Plot the meridional-wind differences.
  (levels, colormap) = symmetric_contours(minimum(deltavls), maximum(deltavls))
  subplot(222)
  meridional_wind = contourf(x, z, deltavls; levels = levels, cmap = colormap)
  xlabel(L"x\,\left[\mathrm{km}\right]")
  ylabel(L"z\,\left[\mathrm{km}\right]")
  title(L"Large-scale meridional wind $v\,\left[\mathrm{m\,s^{-1}}\right]$")
  colorbar(meridional_wind)

  # Plot the vertical-wind differences.
  (levels, colormap) = symmetric_contours(minimum(deltawls), maximum(deltawls))
  subplot(223)
  vertical_wind = contourf(x, z, deltawls; levels = levels, cmap = colormap)
  xlabel(L"x\,\left[\mathrm{km}\right]")
  ylabel(L"z\,\left[\mathrm{km}\right]")
  title(L"Large-scale vertical wind $w\,\left[\mathrm{m\,s^{-1}}\right]$")
  colorbar(vertical_wind)

  # Plot the tracer differences.
  (levels, colormap) = symmetric_contours(minimum(deltatls), maximum(deltatls))
  subplot(224)
  tracer = contourf(x, z, deltatls; levels = levels, cmap = colormap)
  xlabel(L"x\,\left[\mathrm{km}\right]")
  ylabel(L"z\,\left[\mathrm{km}\right]")
  title("Large-scale tracer")
  colorbar(tracer)

  # Save the figure.
  savefig(data_path * "results/2DWP90_msgwam_difference.png")
end
