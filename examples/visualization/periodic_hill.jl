# examples/visualization/periodic_hill.jl

using HDF5
using LaTeXStrings

include("style.jl")

# Import the data.
@ivy if length(ARGS) == 0
    data = h5open("./pincflow_output.h5")
elseif length(ARGS) == 1
    data = h5open(ARGS[1] * "/pincflow_output.h5")
else
    error("Too many arguments to the script!")
end

# Set the grid.
x = @. data["x"][:] / 1000
z = @. data["z"][:, 1, :] / 1000
x = @. x * ones(size(z))

# Get the vertical wind.
w = data["w"][:, 1, :, end]

# Close the file.
close(data)

# Create the plot.
(levels, colormap) = symmetric_contours(minimum(w), maximum(w))
contours = contourf(x, z, w; levels = levels, cmap = colormap)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")
savefig("examples/results/periodic_hill.png")
clf()
