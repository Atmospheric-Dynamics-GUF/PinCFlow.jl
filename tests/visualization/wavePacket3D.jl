include("style.jl")
include("tools.jl")

using LaTeXStrings

# Set paths.
data_path = "../"
reference_path = "../"

println("data_path =", data_path)
println("reference_path =", reference_path)

# Import data.
data = ModelOutput(data_path * "wavePacket3D")
reference = ModelOutput(reference_path * "wavePacket3D")

println("Time: ", data.variables["t"][end], " s")

# Set grid.
x = data.variables["x"][120:390] .* 0.001
z = data.variables["z"][120:420] .* 0.001
x = [xi for xi in x, j in 1:size(z)[1]]
z = [zj for i in 1:size(x)[1], zj in z]

# Get density fluctuations.
rhop = data.groups["atmvar"].variables["rhop"][120:390, 1, 120:420, end]
rhopr = reference.groups["atmvar"].variables["rhop"][120:390, 1, 120:420, end]
deltarhop = rhop .- rhopr

# Create plot.
(levels, colormap) = symmetric_contours(minimum(rhop), maximum(rhop))
contours = contourf(x, z, rhop; levels = levels, cmap = colormap)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
colorbar(contours; label = L"\rho'\,\left[\mathrm{kg\,m^{- 3}}\right]")
savefig(data_path * "/results/wavePacket3D.png")

if data_path != reference_path
  (levels, colormap) =
    symmetric_contours(minimum(deltarhop), maximum(deltarhop))
  clf()
  contours = contourf(x, z, deltarhop; levels = levels, cmap = colormap)
  xlabel(L"x\,\left[\mathrm{km}\right]")
  ylabel(L"z\,\left[\mathrm{km}\right]")
  colorbar(contours; label = L"\Delta\rho'\,\left[\mathrm{kg\,m^{- 3}}\right]")
  savefig(data_path * "/results/wavePacket3D_difference.png")
end
