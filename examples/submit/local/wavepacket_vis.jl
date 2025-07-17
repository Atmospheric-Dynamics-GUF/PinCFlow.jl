include("style.jl")

using HDF5
using LaTeXStrings
using Statistics

# Set paths.
data_path = "./examples/submit/local"

println("data_path =", data_path)

# Import data.
data = h5open(data_path * "/pincflow_output.h5")

println(keys(data))

println("Time: ", data["t"][end], " s")

kzmax = 30

# Set grid.
x = data["x"][:] .* 0.001
y = data["y"][:] .* 0.001
z = data["z"][:, :, 1:Int(kzmax)] .* 0.001
t = data["t"][:]
x = [xi for xi in x, iy in 1:size(z)[2], iz in 1:size(z)[3]]
y = [yi for ix in 1:size(z)[1], yi in y, iz in 1:size(z)[3]]

u = data["u"][:, 1, 1:Int(kzmax), end]
w = data["w"][:, 1, 1:Int(kzmax), end]
chi = data["chi"][:, 1, 1:Int(kzmax), end] - 
    data["chi"][:, 1, 1:Int(kzmax), begin]
dchidt = - data["dchidt"][:, 1, 1:Int(kzmax), end]
uchi = data["uchi"][:, 1, 1:Int(kzmax), end]

# Close files.
close(data)

println("Maximum zonal flux = ", maximum(uchi))
println("Maximum Q0 = ", maximum(abs.(dchidt)))
println("Maximum tracer change = ", maximum(abs.(chi)))

# Create the figure.
figure(; figsize = (16, 5))

subplot(131)
variable = uchi
(levels, colormap) =
    symmetric_contours(-maximum(variable), maximum(variable); number = 40)

contours =
    contourf(x[:, 1, :], z[:, 1, :], variable; levels = levels, cmap = colormap)
xlabel("x [km]")
ylabel("altitude [km]")
title("Zonal flux [m/s]")
colorbar(contours)

subplot(132)
variable = dchidt

(levels, colormap) =
    symmetric_contours(-maximum(variable), maximum(variable); number = 40)
contours =
    contourf(x[:, 1, :], z[:, 1, :], variable; levels = levels, cmap = colormap)
xlabel("x [km]")
ylabel("altitude [km]")
title(L"$Q^{(0)}$")
colorbar(contours)

subplot(133)
variable = chi

(levels, colormap) =
    symmetric_contours(-maximum(variable), maximum(variable); number = 40)
contours =
    contourf(x[:, 1, :], z[:, 1, :], variable; levels = levels, cmap = colormap)
xlabel("x [km]")
ylabel("altitude [km]")
title(L"Change in tracer $\Delta\psi$")
colorbar(contours)

savefig("./examples/submit/local/wavepacket_2DWP90.pdf")