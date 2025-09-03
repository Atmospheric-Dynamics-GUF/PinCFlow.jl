include("style.jl")

using HDF5
using LaTeXStrings

data_path = "./examples/submit/local"

println("data_path =", data_path)

# Import data.
data = h5open(data_path * "/pincflow_output_mountain.h5")

println("Time: ", data["t"][end], " s")

# Set grid.
x = data["x"][:] .* 0.001 .- 10
y = data["y"][:] .* 0.001 .- 10
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Get vertical wind.
# w = data["chi"][:, :, :, end] .- data["chi"][:, :, :, begin]
w = data["w"][:, :, :, end]

# Close files.
close(data)

# Create the figure.
figure(; figsize = (12, 3))

# Plot in x-y plane.
iz = 10
subplot(131)
(levels, colormap) =
    symmetric_contours(minimum(w[:, :, iz]), maximum(w[:, :, iz]))
contours = contourf(
    x[:, :, iz],
    y[:, :, iz],
    w[:, :, iz];
    levels = levels,
    cmap = colormap,
)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"y\,\left[\mathrm{km}\right]")
title(L"z\approx 5\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Plot in x-z plane.
iy = 20
subplot(132)
(levels, colormap) =
    symmetric_contours(minimum(w[:, iy, :]), maximum(w[:, iy, :]))
contours = contourf(
    x[:, iy, :],
    z[:, iy, :],
    w[:, iy, :];
    levels = levels,
    cmap = colormap,
)
plot(x[:, iy, 1], z[:, iy, 1]; color = "black", linewidth = 0.5)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"y\approx 0\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Plot in y-z plane.
ix = 20
subplot(133)
(levels, colormap) =
    symmetric_contours(minimum(w[ix, :, :]), maximum(w[ix, :, :]))
contours = contourf(
    y[ix, :, :],
    z[ix, :, :],
    w[ix, :, :];
    levels = levels,
    cmap = colormap,
)
plot(y[ix, :, 1], z[ix, :, 1]; color = "black", linewidth = 0.5)
xlabel(L"y\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"x\approx 0\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Save the figure.
savefig("./examples/submit/local/mountain_wave.png")
