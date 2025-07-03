include("style.jl")
include("tools.jl")

using LaTeXStrings

# Set paths.
data_path = "../"
reference_path = "../"

println("data_path =", data_path)
println("reference_path =", reference_path)

# Import data.
data = ModelOutput(data_path * "mountainwave")
reference = ModelOutput(reference_path * "mountainwave")

println("Time: ", data.variables["t"][end], " s")

# Set grid.
x = data.variables["x"][99:202] .* 0.001 .- 30
z = data.variables["z"][99:202, 1, 1:51, end] .* 0.001
x = x .* ones(size(z))

# Get vertical wind.
w = data.groups["atmvar"].variables["w"][99:202, 1, 1:51, end]
wr = reference.groups["atmvar"].variables["w"][99:202, 1, 1:51, end]
deltaw = w .- wr

# Create plot.
(levels, colormap) = symmetric_contours(minimum(w), maximum(w))
contours = contourf(x, z, w; levels = levels, cmap = colormap)
plot(x[:, 1], z[:, 1]; color = "black", linewidth = 1.0)
xlim(-10.0, 10.0)
ylim(0.0, 10.0)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")
savefig(data_path * "/results/mountainwave.png")

# Create difference plot.
if data_path != reference_path
  (levels, colormap) = symmetric_contours(minimum(deltaw), maximum(deltaw))
  clf()
  contours = contourf(x, z, deltaw; levels = levels, cmap = colormap)
  plot(x[:, 1], z[:, 1]; color = "black", linewidth = 1.0)
  xlim(-10.0, 10.0)
  ylim(0.0, 10.0)
  xlabel(L"x\,\left[\mathrm{km}\right]")
  ylabel(L"z\,\left[\mathrm{km}\right]")
  colorbar(contours; label = L"\Delta w\,\left[\mathrm{m\,s^{-1}}\right]")
  savefig(data_path * "/results/mountainwave_difference.png")
end
