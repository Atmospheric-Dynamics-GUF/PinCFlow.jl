# examples/visualization/wkb_mountain_wave.jl

using HDF5
using CairoMakie
using LaTeXStrings
using PinCFlow

set_visualization_theme!()

# Import the data.
@ivy if length(ARGS) == 0
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
x = [xi for xi in x, j in 1:size(z)[2], k in 1:size(z)[3]]
y = [yj for i in 1:size(z)[1], yj in y, k in 1:size(z)[3]]

# Get the vertical wind.
w = data["w"][:, :, :, end]

# Close the file.
close(data)

# Create the figure.
figure = Figure()

# Plot in x-y plane.
k = 10
axis = Axis(
    figure[1, 1];
    title = L"z\approx 5\,\mathrm{km}",
    xlabel = L"x\,[\mathrm{km}]",
    ylabel = L"y\,[\mathrm{km}]",
)
@ivy (levels, colormap) =
    symmetric_contours(minimum(w[:, :, k]), maximum(w[:, :, k]))
@ivy contours =
    contourf!(axis, x[:, :, k], y[:, :, k], w[:, :, k]; levels, colormap)
tightlimits!(axis)
Colorbar(
    figure[1, 2],
    contours;
    ticks = trunc.(levels; digits = 4),
    label = L"w\,[\mathrm{m\,s^{-1}}]",
)

# Plot in x-z plane.
j = 20
axis = Axis(
    figure[1, 3];
    title = L"y\approx 0\,\mathrm{km}",
    xlabel = L"x\,[\mathrm{km}]",
    ylabel = L"z\,[\mathrm{km}]",
)
@ivy (levels, colormap) =
    symmetric_contours(minimum(w[:, j, :]), maximum(w[:, j, :]))
@ivy contours =
    contourf!(axis, x[:, j, :], z[:, j, :], w[:, j, :]; levels, colormap)
tightlimits!(axis)
@ivy lines(x[:, j, 1], z[:, j, 1]; color = :black, linewidth = 0.5)
Colorbar(
    figure[1, 4],
    contours;
    ticks = trunc.(levels; digits = 4),
    label = L"w\,[\mathrm{m\,s^{-1}}]",
)

# Plot in y-z plane.
i = 20
axis = Axis(
    figure[1, 5];
    title = L"x\approx 0\,\mathrm{km}",
    xlabel = L"y\,[\mathrm{km}]",
    ylabel = L"z\,[\mathrm{km}]",
)
@ivy (levels, colormap) =
    symmetric_contours(minimum(w[i, :, :]), maximum(w[i, :, :]))
@ivy contours =
    contourf!(axis, y[i, :, :], z[i, :, :], w[i, :, :]; levels, colormap)
tightlimits!(axis)
@ivy lines(y[i, :, 1], z[i, :, 1]; color = :black, linewidth = 0.5)
Colorbar(
    figure[1, 6],
    contours;
    ticks = trunc.(levels; digits = 4),
    label = L"w\,[\mathrm{m\,s^{-1}}]",
)

# Resize, display and save the figure.
resize_to_layout!(figure)
display(figure)
save("examples/results/wkb_mountain_wave.svg", figure)
