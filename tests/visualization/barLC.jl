include("style.jl")
include("tools.jl")

using FFTW
using LaTeXStrings

# Set paths.
data_path = "../"
reference_path = "../"

println("data_path =", data_path)
println("reference_path =", reference_path)

# Import data.
data = ModelOutput(data_path * "barLC")
reference = ModelOutput(reference_path * "barLC")

println("Time: ", data.variables["t"][end], " s")

# Set grid.
x = data.variables["x"][:] .* 0.001
y = data.variables["y"][43:end] .* 0.001
dx = abs(x[2] - x[1]) * 1000.0
dy = abs(y[2] - y[1]) * 1000.0
x = [xi for xi in x, j in 1:size(y)[1]]
y = [yj for i in 1:size(x)[1], yj in y]

# Loop over data and reference.
ugravity = Dict()
thetap = Dict()
for (path, data_set) in zip((reference_path, data_path), (reference, data))

  # Set fields of interest.
  u = data_set.groups["atmvar"].variables["u"][:, 43:end, 3, end]
  v = data_set.groups["atmvar"].variables["v"][:, 43:end, 3, end]
  thetap[path] =
    data_set.groups["atmvar"].variables["thetap"][:, 43:end, 3, end]

  # Compute divergence.
  divergence = zeros(size(u))
  for j in 2:(size(u)[2] - 1), i in 2:(size(u)[1] - 1)
    divergence[i, j] =
      0.5 *
      ((u[i + 1, j] - u[i - 1, j]) / dx + (v[i, j + 1] - v[i, j - 1]) / dy)
  end

  # Apply Fourier filter.
  sigma = divergence
  sigmatilde = fft(sigma)
  k = fftfreq(size(sigma)[1], 1.0 / dx)
  l = fftfreq(size(sigma)[2], 1.0 / dy)
  urossbytilde = deepcopy(sigmatilde)
  urossbytilde[:, abs.(l) .> 1.0e-6] .= 0.0
  urossbytilde[abs.(k) .> 1.0e-6, :] .= 0.0
  urossby = real.(ifft(urossbytilde))
  ugravity[path] = sigma .- urossby
end

# Compute difference of relevant fields.
deltaugravity = ugravity[data_path] .- ugravity[reference_path]
deltathetap = thetap[data_path] .- thetap[reference_path]
ugravity = ugravity[data_path]
thetap = thetap[data_path]

# Create plot.
(levels, colormap) = symmetric_contours(minimum(ugravity), maximum(ugravity))
contours = contourf(x, y, ugravity; levels = levels, cmap = colormap)
contour(x, y, thetap; linewidths = 1.0, colors = "black")
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"y\,\left[\mathrm{km}\right]")
colorbar(
  contours;
  label = L"\boldsymbol{\nabla}_z\cdot\boldsymbol{u}\," *
          L"\left[\mathrm{s^{-1}}\right]",
)
savefig(data_path * "/results/barLC.png")

# Create difference plot.
if data_path != reference_path
  (levels, colormap) =
    symmetric_contours(minimum(deltaugravity), maximum(deltaugravity))
  clf()
  contours = contourf(x, y, deltaugravity; levels = levels, cmap = colormap)
  contour(x, y, deltathetap; linewidths = 1.0, colors = "black")
  xlabel(L"x\,\left[\mathrm{km}\right]")
  ylabel(L"y\,\left[\mathrm{km}\right]")
  colorbar(
    contours;
    label = L"\Delta\boldsymbol{\nabla}_z\cdot\boldsymbol{u}" *
            L"\,\left[\mathrm{s^{-1}}\right]",
  )
  savefig(data_path * "/results/barLC_difference.png")
end
