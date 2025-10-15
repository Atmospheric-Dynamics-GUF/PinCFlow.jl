# Mountain-wave simulation

## Simulation

The script

```julia
# examples/submit/mountain_wave.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

h0 = 100.0
l0 = 1000.0

lx = 20000.0
ly = 20000.0
lz = 20000.0
dxr = lx / 2
dyr = ly / 2
dzr = lz / 2
alpharmax = 0.0179

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
)
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 40,
    lx,
    ly,
    lz,
    npx,
    npy,
    npz,
)
grid = GridNamelist(;
    resolved_topography = (x, y) -> h0 / (1 + (x^2 + y^2) / l0^2),
)
output =
    OutputNamelist(; output_variables = (:w,), output_file = "mountain_wave.h5")
sponge = SpongeNamelist(;
    lhs_sponge = (x, y, z, t, dt) -> begin
        alpharx =
            abs(x) >= (lx - dxr) / 2 ?
            sin(pi * (abs(x) - (lx - dxr) / 2) / dxr)^2 : 0.0
        alphary =
            abs(y) >= (ly - dyr) / 2 ?
            sin(pi * (abs(y) - (ly - dyr) / 2) / dyr)^2 : 0.0
        alpharz =
            z >= lz - dzr ? sin(pi / 2 * (z - (lz - dzr)) / dzr)^2 : 0.0
        return alpharmax * (alpharx + alphary + alpharz) / 3
    end,
    relaxed_u = (x, y, z, t, dt) -> 10.0,
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))

```

performs a 3D mountain-wave simulation in 64 MPI processes (4 for each dimension of physical space) and writes the vertical wind to `mountain_wave.h5`, if executed with

```shell
mpiexec=$(julia --project=examples -e 'using MPICH_jll; println(MPICH_jll.mpiexec_path)')
${mpiexec} -n 64 julia examples/submit/mountain_wave.jl 4 4 4
```

(provided the examples project has been set up as illustrated in the user guide and MPI.jl and HDF5.jl are configured to use their default backends). The surface topography is given by

$$h \left(x, y\right) = \frac{h_0}{1 + \left(x^2 + y^2\right) / l_0^2},$$

with $h_0 = 100 \, \mathrm{m}$ and $l_0 = 1 \, \mathrm{km}$. The atmosphere is isothermal, with the default temperature $T_0 = 300 \, \mathrm{K}$ and the initial wind $\boldsymbol{u}_0 = \left(10, 0, 0\right)^\mathrm{T} \, \mathrm{m \, s^{- 1}}$.

Reflections at the upper boundary are prevented by damping the generated mountain waves in a sponge defined by

$$\alpha_\mathrm{R} \left(x, y, z\right) = \alpha_{\mathrm{R}, \max} \frac{\alpha_{\mathrm{R}, x} \left(x\right) + \alpha_{\mathrm{R}, y} \left(y\right) + \alpha_{\mathrm{R}, z} \left(z\right)}{3}$$

with

$$\begin{align*}
    \alpha_{\mathrm{R}, x} \left(x\right) & = \begin{cases}
        \sin^2 \left[\pi \frac{\left|x\right| - \left(L_x - \Delta x_\mathrm{R}\right) / 2}{\Delta x_\mathrm{R}}\right] & \mathrm{if} \quad \left|x\right| \geq \frac{1}{2} \left(L_x - \Delta x_\mathrm{R}\right),\\
        0 & \mathrm{else},
    \end{cases}\\
    \alpha_{\mathrm{R}, y} \left(y\right) & = \begin{cases}
        \sin^2 \left[\pi \frac{\left|y\right| - \left(L_y - \Delta y_\mathrm{R}\right) / 2}{\Delta y_\mathrm{R}}\right] & \mathrm{if} \quad \left|y\right| \geq \frac{1}{2} \left(L_y - \Delta y_\mathrm{R}\right),\\
        0 & \mathrm{else},
    \end{cases}\\
    \alpha_{\mathrm{R}, z} \left(z\right) & = \begin{cases}
        \sin^2 \left[\frac{\pi}{2} \frac{z - \left(L_z - \Delta z_\mathrm{R}\right)}{\Delta z_\mathrm{R}}\right] & \mathrm{if} \quad z \geq L_z - \Delta z_\mathrm{R},\\
        0 & \mathrm{else},
    \end{cases}
\end{align*}$$

where $\alpha_{\mathrm{R}, \max} = 0.0179 \, \mathrm{s^{- 1}}$, $\Delta x_\mathrm{R} = L_x / 2$, $\Delta y_\mathrm{R} = L_y / 2$ and $\Delta z_\mathrm{R} = L_z / 2$. This sponge not only prevents wave reflections at the model top but also provides a damping at the horizontal boundaries. Moreover, it is configured such that the wind is relaxed towards its initial state, so that (in the ideal case) the periodicity in $x$ and $y$ is effectively eliminated by enforcing a constant wind at the domain edges.

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

h5open("mountain_wave.h5") do data
    plot_contours(
        "examples/results/mountain_wave.svg",
        data,
        "w",
        (20, 20, 10, 2);
        label = L"w\,[\mathrm{m\,s^{-1}}]",
    )
    return
end

```

visualizes the vertical wind at the end of the above simulation (i.e. after one hour) in three cross sections of the domain and saves the generated figure to an SVG file that is included below.

![](results/mountain_wave.svg)
