include("style.jl")
include("tools.jl")

using LaTeXStrings

# Set paths.
data_path = "../"
reference_path = "../"

println("data_path =", data_path)
println("reference_path =", reference_path)

# Import data.
data = ModelOutput(data_path * "IGW")
reference = ModelOutput(reference_path * "IGW")

println("Time: ", data.variables["t"][end], " s")

# Set grid.
x = data.variables["x"][:] .* 0.001
z = data.variables["z"][:] .* 0.001
x = [xi for xi in x, j in 1:size(z)[1]]
z = [zj for i in 1:size(x)[1], zj in z]

# Get potential-temperature fluctuations.
thetap = data.groups["atmvar"].variables["thetap"][:, 1, :, end]
thetapr = reference.groups["atmvar"].variables["thetap"][:, 1, :, end]
deltathetap = thetap .- thetapr

# Create plot.
(levels, colormap) = symmetric_contours(minimum(thetap), maximum(thetap))
contours = contourf(x, z, thetap; levels = levels, cmap = colormap)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
colorbar(contours; label = L"\theta'\,\left[\mathrm{K}\right]")
savefig(data_path * "/results/IGW.png")

# Create difference plot.
if data_path != reference_path
  (levels, colormap) =
    symmetric_contours(minimum(deltathetap), maximum(deltathetap))
  clf()
  contours = contourf(x, z, deltathetap; levels = levels, cmap = colormap)
  xlabel(L"x\,\left[\mathrm{km}\right]")
  ylabel(L"z\,\left[\mathrm{km}\right]")
  colorbar(contours; label = L"\Delta\theta'\,\left[\mathrm{K}\right]")
  savefig(data_path * "/results/IGW_difference.png")
end
