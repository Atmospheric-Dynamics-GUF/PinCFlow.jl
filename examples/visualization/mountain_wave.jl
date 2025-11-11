# examples/visualization/mountain_wave.jl

using HDF5
using PythonPlot
using LaTeXStrings
using PinCFlow

<<<<<<< HEAD
# Set paths.
host_name = readchomp(`hostname`)
user_name = readchomp(`whoami`)
if occursin("login", host_name)
	data_path =
		"/scratch/atmodynamics/" * user_name * "/pinc/examples/mountain_wave/"
	reference_path = data_path
elseif occursin("dkrz", host_name)
	data_path = "/scratch/b/" * user_name * "/pinc/examples/mountain_wave/"
	reference_path = data_path
=======
set_plot_style()

# Import the data.
@ivy if length(ARGS) == 0
    data = h5open("./pincflow_output.h5")
elseif length(ARGS) == 1
    data = h5open(ARGS[1] * "/pincflow_output.h5")
>>>>>>> 2aee3f7
else
	data_path = "~/PF/pinc/test/"
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
figure(; figsize = (12, 3))

# Plot in x-y plane.
k = 10
subplot(131)
<<<<<<< HEAD
(levels, colormap) =
	symmetric_contours(minimum(w[:, :, iz]), maximum(w[:, :, iz]))
contours = contourf(
	x[:, :, iz],
	y[:, :, iz],
	w[:, :, iz];
	levels = levels,
	cmap = colormap,
=======
@ivy (levels, colormap) =
    symmetric_contours(minimum(w[:, :, k]), maximum(w[:, :, k]))
@ivy contours = contourf(
    x[:, :, k],
    y[:, :, k],
    w[:, :, k];
    levels = levels,
    cmap = colormap,
>>>>>>> 2aee3f7
)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"y\,\left[\mathrm{km}\right]")
title(L"z\approx 5\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Plot in x-z plane.
j = 20
subplot(132)
<<<<<<< HEAD
(levels, colormap) =
	symmetric_contours(minimum(w[:, iy, :]), maximum(w[:, iy, :]))
contours = contourf(
	x[:, iy, :],
	z[:, iy, :],
	w[:, iy, :];
	levels = levels,
	cmap = colormap,
=======
@ivy (levels, colormap) =
    symmetric_contours(minimum(w[:, j, :]), maximum(w[:, j, :]))
@ivy contours = contourf(
    x[:, j, :],
    z[:, j, :],
    w[:, j, :];
    levels = levels,
    cmap = colormap,
>>>>>>> 2aee3f7
)
@ivy plot(x[:, j, 1], z[:, j, 1]; color = "black", linewidth = 0.5)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"y\approx 0\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Plot in y-z plane.
i = 20
subplot(133)
<<<<<<< HEAD
(levels, colormap) =
	symmetric_contours(minimum(w[ix, :, :]), maximum(w[ix, :, :]))
contours = contourf(
	y[ix, :, :],
	z[ix, :, :],
	w[ix, :, :];
	levels = levels,
	cmap = colormap,
=======
@ivy (levels, colormap) =
    symmetric_contours(minimum(w[i, :, :]), maximum(w[i, :, :]))
@ivy contours = contourf(
    y[i, :, :],
    z[i, :, :],
    w[i, :, :];
    levels = levels,
    cmap = colormap,
>>>>>>> 2aee3f7
)
@ivy plot(y[i, :, 1], z[i, :, 1]; color = "black", linewidth = 0.5)
xlabel(L"y\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
title(L"x\approx 0\,\mathrm{km}")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")

# Save the figure.
savefig("examples/results/mountain_wave.png")
clf()
