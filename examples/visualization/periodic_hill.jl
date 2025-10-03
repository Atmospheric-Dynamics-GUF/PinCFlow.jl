# examples/visualization/periodic_hill.jl

using HDF5
using CairoMakie
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
z = data["z"][:, 1, :] ./ 1000
x = x .* ones(size(z))

# Get the vertical wind.
w = data["w"][:, 1, :, end]

# Close the file.
close(data)

# Create the plot.
(levels, colormap) = symmetric_contours(minimum(w), maximum(w))
figure = Figure()
axis = Axis(
    figure[1, 1];
    xlabel = L"x\,[\mathrm{km}]",
    ylabel = L"z\,[\mathrm{km}]",
)
contours = contourf!(axis, x, z, w; levels, colormap)
tightlimits!(axis)
Colorbar(
    figure[1, 2],
    contours;
    ticks = trunc.(levels; digits = 4),
    label = L"w\,[\mathrm{m\,s^{-1}}]",
)
resize_to_layout!(figure)
display(figure)
save("examples/results/periodic_hill.svg", figure)
