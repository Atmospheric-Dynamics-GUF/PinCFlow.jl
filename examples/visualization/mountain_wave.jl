# examples/visualization/mountain_wave.jl

using HDF5
using LaTeXStrings

include("style.jl")

# Import the data.
if length(ARGS) == 0
    data = h5open("./pincflow_output.h5")
elseif length(ARGS) == 1
    data = h5open(ARGS[1] * "/pincflow_output.h5")
else
    error("Too many arguments to the script!")
end

# Set the grid.
x = data["x"][:] ./ 1000
y = data["y"][:] ./ 1000
z = data["z"][:, :, :] ./ 1000
x = [xi for xi in x, j in size(z)[2], k in size(z)[3]]
y = [yj for i in size(z)[1], yj in y, k in size(z)[3]]

# Get the vertical wind.
w = data["w"][:, :, :, end]

# Close the file.
close(data)

# Create the figure.
figure(; figsize = (12, 3))

# Plot in x-y plane.
k = 10
subplot(131)
@views (levels, colormap) =
    symmetric_contours(minimum(w[:, :, k]), maximum(w[:, :, k]))
@views contours = contourf(
    x[:, :, k],
    y[:, :, k],
    w[:, :, k];
    levels = levels,
    cmap = colormap,
)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"y\,\left[\mathrm{km}\right]")
title(L"z\approx 5\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Plot in x-z plane.
j = 20
subplot(132)
@views (levels, colormap) =
    symmetric_contours(minimum(w[:, j, :]), maximum(w[:, j, :]))
@views contours = contourf(
    x[:, j, :],
    z[:, j, :],
    w[:, j, :];
    levels = levels,
    cmap = colormap,
)
@views plot(x[:, j, 1], z[:, j, 1]; color = "black", linewidth = 0.5)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"y\approx 0\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Plot in y-z plane.
i = 20
subplot(133)
@views (levels, colormap) =
    symmetric_contours(minimum(w[i, :, :]), maximum(w[i, :, :]))
@views contours = contourf(
    y[i, :, :],
    z[i, :, :],
    w[i, :, :];
    levels = levels,
    cmap = colormap,
)
@views plot(y[i, :, 1], z[i, :, 1]; color = "black", linewidth = 0.5)
xlabel(L"y\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"x\approx 0\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Save the figure.
savefig("examples/results/mountain_wave.png")
clf()
