include("style.jl")

using HDF5
using LaTeXStrings

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
end

println("data_path =", data_path)
println("reference_path =", reference_path)

# Import data.
data = h5open(data_path * "/pincflow_output.h5")
reference = h5open(reference_path * "/pincflow_output.h5")

println("Time: ", data["t"][end], " s")

# Set grid.
x = data["x"][:] .* 0.001 .- 10
y = data["y"][:] .* 0.001 .- 10
z = data["z"][:, :, :] .* 0.001
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

# Get vertical wind.
w = data["w"][:, :, :, end]
wr = reference["w"][:, :, :, end]
deltaw = w .- wr

# Close files.
close(data)
close(reference)

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
savefig("../results/mountain_wave.png")

# Create difference plots.
if reference_path != data_path

    # Create the figure.
    figure(; figsize = (12, 3))

    # Plot in x-y plane.
    iz = 10
    subplot(131)
    (levels, colormap) =
        symmetric_contours(minimum(deltaw[:, :, iz]), maximum(deltaw[:, :, iz]))
    contours = contourf(
        x[:, :, iz],
        y[:, :, iz],
        deltaw[:, :, iz];
        levels = levels,
        cmap = colormap,
    )
    xlabel(L"x\,\left[\mathrm{km}\right]")
    ylabel(L"y\,\left[\mathrm{km}\right]")
    title(L"z\approx 5\,\mathrm{km}")
    colorbar(contours; label = L"\Delta w\,\left[\mathrm{m\,s^{-1}}\right]")

    # Plot in x-z plane.
    iy = 20
    subplot(132)
    (levels, colormap) =
        symmetric_contours(minimum(deltaw[:, iy, :]), maximum(deltaw[:, iy, :]))
    contours = contourf(
        x[:, iy, :],
        z[:, iy, :],
        deltaw[:, iy, :];
        levels = levels,
        cmap = colormap,
    )
    plot(x[:, iy, 1], z[:, iy, 1]; color = "black", linewidth = 0.5)
    xlabel(L"x\,\left[\mathrm{km}\right]")
    ylabel(L"z\,\left[\mathrm{km}\right]")
    title(L"y\approx 0\,\mathrm{km}")
    colorbar(contours; label = L"\Delta w\,\left[\mathrm{m\,s^{-1}}\right]")

    # Plot in y-z plane.
    ix = 20
    subplot(133)
    (levels, colormap) =
        symmetric_contours(minimum(deltaw[ix, :, :]), maximum(deltaw[ix, :, :]))
    contours = contourf(
        y[ix, :, :],
        z[ix, :, :],
        deltaw[ix, :, :];
        levels = levels,
        cmap = colormap,
    )
    plot(y[ix, :, 1], z[ix, :, 1]; color = "black", linewidth = 0.5)
    xlabel(L"y\,\left[\mathrm{km}\right]")
    ylabel(L"z\,\left[\mathrm{km}\right]")
    title(L"x\approx 0\,\mathrm{km}")
    colorbar(contours; label = L"\Delta w\,\left[\mathrm{m\,s^{-1}}\right]")

    # Save the figure.
    savefig("../results/mountain_wave_differences.png")
end
