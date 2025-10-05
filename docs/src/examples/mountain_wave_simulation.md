# Mountain-wave simulation

## Simulation

The script

```julia
# examples/submit/mountain_wave.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

atmosphere = AtmosphereNamelist(; initial_wind = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 40,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
    npx = 8,
    npy = 8,
)
grid = GridNamelist(; mountain_case = 4)
output =
    OutputNamelist(; output_variables = (:w,), output_file = "mountain_wave.h5")
sponge = SpongeNamelist(;
    use_sponge = true,
    alpharmax = 1.79E-2,
    betarmax = 0.0E+0,
    lateral_sponge = true,
    sponge_type = SinusoidalSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))

```

performs a 3D mountain-wave simulation with parallelization in the zonal and meridional dimensions, and writes the vertical wind to `mountain_wave.h5`, if executed with

```shell
mpiexec=$(julia --project=examples -e 'using MPICH_jll; println(MPICH_jll.mpiexec_path)')
${mpiexec} -n 64 julia --project examples/submit/mountain_wave.jl
```

from the root directory of the repository (provided MPI.jl and HDF5.jl are configured to use their default backends). The surface topography is given by

$$h \left(x, y\right) = \frac{h_0}{1 + \left(x^2 + y^2\right) / l_0^2},$$

where the default values $h_0 = 100 \, \mathrm{m}$ and $l_0 = 1 \, \mathrm{km}$ are being used. The atmosphere is isothermal, with the default temperature $T_0 = 300 \, \mathrm{K}$ and the initial wind $\boldsymbol{u}_0 = \left(10, 0, 0\right)^\mathrm{T} \, \mathrm{m \, s^{- 1}}$.

Reflections at the upper boundary are prevented by damping the generated mountain waves in a sponge defined by

$$\alpha_\mathrm{R} \left(x, y, z\right) = \frac{\alpha_{\mathrm{R}, x} \left(x\right) + \alpha_{\mathrm{R}, y} \left(y\right) + \alpha_{\mathrm{R}, z} \left(z\right)}{3}$$

with

$$\begin{align*}
    \alpha_{\mathrm{R}, x} \left(x\right) & = \begin{cases}
        \alpha_{\mathrm{R}, \max} \sin^2 \left[\frac{\pi \left(x_{\mathrm{R}, 0} - x\right)}{2 \Delta x_\mathrm{R}}\right] & \mathrm{if} \quad x \leq x_{\mathrm{R}, 0},\\
        \alpha_{\mathrm{R}, \max} \sin^2 \left[\frac{\pi \left(x - x_{\mathrm{R}, 1}\right)}{2 \Delta x_\mathrm{R}}\right] & \mathrm{if} \quad x \geq x_{\mathrm{R}, 1},\\
        0 & \mathrm{else},
    \end{cases}\\
    \alpha_{\mathrm{R}, y} \left(y\right) & = \begin{cases}
        \alpha_{\mathrm{R}, \max} \sin^2 \left[\frac{\pi \left(y_{\mathrm{R}, 0} - y\right)}{2 \Delta y_\mathrm{R}}\right] & \mathrm{if} \quad y \leq y_{\mathrm{R}, 0},\\
        \alpha_{\mathrm{R}, \max} \sin^2 \left[\frac{\pi \left(y - y_{\mathrm{R}, 1}\right)}{2 \Delta y_\mathrm{R}}\right] & \mathrm{if} \quad y \geq y_{\mathrm{R}, 1},\\
        0 & \mathrm{else},
    \end{cases}\\
    \alpha_{\mathrm{R}, z} \left(z\right) & = \begin{cases}
        \alpha_{\mathrm{R}, \max} \sin^2 \left[\frac{\pi \left(z - z_\mathrm{R}\right)}{2 \Delta z_\mathrm{R}}\right] & \mathrm{if} \quad z \geq z_\mathrm{R},\\
        0 & \mathrm{else},
    \end{cases}
\end{align*}$$

where the maximum of the damping coefficient is $\alpha_{\mathrm{R}, \max} = 0.0179 \, \mathrm{s^{- 1}}$, which corresponds to the buoyancy frequency. Since the simulation uses the default setting `sponge_extent = 0.5`, the parameter $\Delta z_\mathrm{R}$ is given by half of the domain's vertical extent, whereas $\Delta x_\mathrm{R}$ and $\Delta y_\mathrm{R}$ are each given by a quarter of the domain's extent in the respective dimension. The edges of the sponge are such that it is horizontally centered at $\left(- 10, - 10\right)^\mathrm{T} \, \mathrm{km}$ and has an extent of $\left(\Delta x_\mathrm{R}, \Delta y_\mathrm{R}\right)^\mathrm{T}$ below $z_\mathrm{R} = 10 \, \mathrm{km}$, whereas it covers the entire horizontal plane above that altitude (see below for plots of $\alpha_\mathrm{R}$ in three cross sections of the domain). This means that the sponge not only prevents wave reflections at the model top but also provides a damping at the horizontal boundaries. Moreover, it is configured such that the wind is relaxed towards its initial state, so that (in the ideal case) the periodicity in $x$ and $y$ is effectively eliminated by enforcing a constant wind at the domain edges.

![](sinusoidal_sponge.svg)

## Visualization

The script

```julia
# examples/visualization/mountain_wave.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

set_visualization_theme!()

# Import the data.
data = h5open("mountain_wave.h5")

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
save("examples/results/mountain_wave.svg", figure)

```

visualizes the vertical wind at the end of the above simulation (i.e. after one hour) in three cross sections of the domain and saves the generated figure to an SVG file that is included below. Note that `symmetric_contours` returns a cropped colormap that is centered at $w = 0 \, \mathrm{m \, s^{- 1}}$.

![](results/mountain_wave.svg)

## See also

  - [`PinCFlow.Types.NamelistTypes.AtmosphereNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.GridNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.SpongeNamelist`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.compute_topography`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Sponge`](@ref)

  - [`PinCFlow.Update.compute_sponges!`](@ref)

  - [`PinCFlow.Update.apply_lhs_sponge!`](@ref)

  - [`PinCFlow.set_visualization_theme!`](@ref)

  - [`PinCFlow.symmetric_contours`](@ref)
