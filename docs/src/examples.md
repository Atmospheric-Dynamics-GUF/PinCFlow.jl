# Examples

## Cold bubble

The script

```julia
# examples/scripts/cold_bubble.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

lx = 20000.0
lz = 20000.0

rx = lx / 8
rz = lz / 8

atmosphere = AtmosphereNamelist(;
    background = Isentropic(),
    initial_rhop = (x, y, z) -> begin
        r = sqrt((x / rx)^2 + ((z - 3 * rz) / rz)^2)
        if r <= 1
            return 0.005 * (1 + cos(pi * r))
        else
            return 0.0
        end
    end,
)
discretization = DiscretizationNamelist(; dtmax = 60.0)
domain = DomainNamelist(; x_size = 40, z_size = 40, lx, lz, npx, npz)
output = OutputNamelist(;
    output_variables = (:thetap,),
    output_file = "cold_bubble.h5",
)

integrate(Namelists(; atmosphere, discretization, domain, output))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("cold_bubble.h5") do data
        plot_output(
            "examples/results/cold_bubble.svg",
            data,
            ("thetap", 1, 1, 1, 2);
        )
        return
    end
end

```

simulates a cold bubble in a 2D pseudo-incompressible isentropic atmosphere and visualizes the potential-temperature fluctuations after one hour integration time (see below).

![](examples/results/cold_bubble.svg)

## Hot bubble

The script

```julia
# examples/scripts/hot_bubble.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

lx = 20000.0
lz = 20000.0

rx = lx / 8
rz = lz / 8

atmosphere = AtmosphereNamelist(;
    model = Compressible(),
    background = Isentropic(),
    initial_rhop = (x, y, z) -> begin
        r = sqrt((x / rx)^2 + ((z - 5 * rz) / rz)^2)
        if r <= 1
            return -0.005 * (1 + cos(pi * r))
        else
            return 0.0
        end
    end,
)
discretization = DiscretizationNamelist(; dtmax = 60.0)
domain = DomainNamelist(; x_size = 40, z_size = 40, lx, lz, npx, npz)
output = OutputNamelist(;
    output_variables = (:thetap,),
    output_file = "hot_bubble.h5",
)

integrate(Namelists(; atmosphere, discretization, domain, output))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("hot_bubble.h5") do data
        plot_output(
            "examples/results/hot_bubble.svg",
            data,
            ("thetap", 1, 1, 1, 2);
        )
        return
    end
end

```

simulates a hot bubble in a 2D compressible isentropic atmosphere and visualizes the potential-temperature fluctuations after one hour integration time (see below).

![](examples/results/hot_bubble.svg)

## Mountain wave

The script

```julia
# examples/scripts/mountain_wave.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
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

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("mountain_wave.h5") do data
        plot_output(
            "examples/results/mountain_wave.svg",
            data,
            ("w", 20, 20, 10, 2);
        )
        return
    end
end

```

performs a 3D mountain-wave simulation. The surface topography is given by

$$h \left(x, y\right) = \frac{h_0}{1 + \left(x^2 + y^2\right) / l_0^2},$$

with $h_0 = 100 \ \mathrm{m}$ and $l_0 = 1 \ \mathrm{km}$. The atmosphere is isothermal, with the default temperature $T_0 = 300 \ \mathrm{K}$ and the initial wind $\boldsymbol{u}_0 = \left(10, 0, 0\right)^\mathrm{T} \ \mathrm{m \ s^{- 1}}$.

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

where $\alpha_{\mathrm{R}, \max} = 0.0179 \ \mathrm{s^{- 1}}$, $\Delta x_\mathrm{R} = L_x / 2$, $\Delta y_\mathrm{R} = L_y / 2$ and $\Delta z_\mathrm{R} = L_z / 2$. This sponge not only prevents wave reflections at the model top but also provides a damping at the horizontal boundaries. Moreover, it is configured such that the wind is relaxed towards its initial state, so that (in the ideal case) the periodicity in $x$ and $y$ is effectively eliminated by enforcing a constant wind at the domain edges.

After the simulation has finished, the vertical wind is visualized in three cross sections of the domain (see below).

![](examples/results/mountain_wave.svg)

## Vortex

The script

```julia
# examples/scripts/vortex.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

lx = 20000.0
ly = 20000.0

rx = lx / 4
ry = ly / 4

atmosphere = AtmosphereNamelist(;
    model = Boussinesq(),
    background = NeutralStratification(),
    initial_u = (x, y, z) -> begin
        r = sqrt((x / rx)^2 + (y / ry)^2)
        if r <= 1
            return -5 * y / ry * (1 + cos(pi * r)) / 2
        else
            return 0.0
        end
    end,
    initial_v = (x, y, z) -> begin
        r = sqrt((x / rx)^2 + (y / ry)^2)
        if r <= 1
            return 5 * x / rx * (1 + cos(pi * r)) / 2
        else
            return 0.0
        end
    end,
)
domain = DomainNamelist(; x_size = 40, y_size = 40, lx, ly, npx, npy)
output = OutputNamelist(;
    output_variables = (:chi, :u, :v),
    output_file = "vortex.h5",
)
tracer = TracerNamelist(;
    tracer_setup = TracerOn(),
    initial_tracer = (x, y, z) -> begin
        r = sqrt(((abs(x) - rx) / rx)^2 + (y / ry)^2)
        if r <= 1
            return sign(x) * (1 + cos(pi * r)) / 2
        else
            return 0.0
        end
    end,
)

integrate(Namelists(; atmosphere, domain, output, tracer))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("vortex.h5") do data
        plot_output("examples/results/vortex.svg", data, ("chi", 1, 1, 1, 2);)
        return
    end
end

```

initializes two tracer disks and a vortex in a 2D horizontal Boussinesq atmosphere, integrates over one hour and visualizes the resulting tracer distribution (see below).

![](examples/results/vortex.svg)

## Wave packet

The script

```julia
# examples/scripts/wave_packet.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

x_size = 40
y_size = 40
z_size = 80

lx = 20000.0
ly = 20000.0
lz = 40000.0

rx = 0.25
ry = 0.25
rz = 0.25

x0 = 0.0
y0 = 0.0
z0 = 20000.0

a0 = 0.05

k = 16 * pi / lx
l = 16 * pi / ly
m = 32 * pi / lz

background = Realistic()
coriolis_frequency = 0.0001

atmosphere = AtmosphereNamelist(; background, coriolis_frequency)
domain = DomainNamelist(;
    x_size,
    y_size,
    z_size,
    lx,
    ly,
    lz,
    base_comm = MPI.COMM_SELF,
)
auxiliary_state = State(Namelists(; atmosphere, domain))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

include("wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(;
    background,
    coriolis_frequency,
    initial_rhop = (x, y, z) ->
        rhobar(x, y, z) *
        (1 / (1 + real(bhat(x, y, z) * exp(1im * phi(x, y, z))) / g) - 1),
    initial_u = (x, y, z) -> real(uhat(x, y, z) * exp(1im * phi(x, y, z))),
    initial_v = (x, y, z) -> real(vhat(x, y, z) * exp(1im * phi(x, y, z))),
    initial_w = (x, y, z) -> real(what(x, y, z) * exp(1im * phi(x, y, z))),
    initial_pip = (x, y, z) ->
        real(pihat(x, y, z) * exp(1im * phi(x, y, z))),
)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :v, :w),
    output_file = "wave_packet.h5",
    tmax = 900.0,
    output_interval = 900.0,
)

integrate(Namelists(; atmosphere, domain, output))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wave_packet.h5") do data
        plot_output(
            "examples/results/wave_packet.svg",
            data,
            ("u", 20, 20, 40, 2),
            ("v", 20, 20, 40, 2),
            ("w", 20, 20, 40, 2);
            time_unit = "min",
        )
        return
    end
end

```

initializes a resolved gravity-wave packet in the stratosphere of a "realistic" atmosphere (isentropic troposphere and isothermal stratosphere) and visualizes the resulting wind after fifteen minutes integration time (see below). For the relatively complex initialization, this script first constructs an auxiliary state that contains the necessary background fields and then uses helper functions that implement the gravity-wave dispersion and polarization relations (included in a separate section below).

![](examples/results/wave_packet.svg)

## WKB mountain wave

The script

```julia
# examples/scripts/wkb_mountain_wave.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

h0 = 150.0
l0 = 5000.0
rl = 10
rh = 2

lx = 400000.0
ly = 400000.0
lz = 20000.0
dxr = lx / 20
dyr = ly / 20
dzr = lz / 10
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
    resolved_topography = (x, y) ->
        x^2 + y^2 <= (rl * l0)^2 ?
        h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) * rh / (rh + 1) : 0.0,
    unresolved_topography = (alpha, x, y) ->
        x^2 + y^2 <= (rl * l0)^2 ?
        (
            pi / l0,
            0.0,
            h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) / (rh + 1),
        ) : (0.0, 0.0, 0.0),
)
output = OutputNamelist(;
    save_ray_volumes = true,
    output_file = "wkb_mountain_wave.h5",
)
sponge = SpongeNamelist(;
    lhs_sponge = (x, y, z, t, dt) ->
        alpharmax / 3 * (
            exp((abs(x) - lx / 2) / dxr) +
            exp((abs(y) - ly / 2) / dyr) +
            exp((z - lz) / dzr)
        ),
    relaxed_u = (x, y, z, t, dt) -> 10.0,
)
wkb = WKBNamelist(; wkb_mode = MultiColumn())

integrate(Namelists(; atmosphere, domain, grid, output, sponge, wkb))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wkb_mountain_wave.h5") do data
        plot_output(
            "examples/results/wkb_mountain_wave.svg",
            data,
            ("nr", 20, 20, 10, 2);
        )
        return
    end
end

```

performs a 3D WKB mountain-wave simulation. The full surface topography is given by

$$\begin{align*}
    h \left(x, y\right) & = \begin{cases}
        \frac{h_0}{2 \left(r_h + 1\right)} \left[1 + \cos \left(\frac{\pi}{r_l l_0} \sqrt{x^2 + y^2}\right)\right] \left[r_h + \cos \left(\frac{\pi x}{l_0}\right)\right] & \mathrm{if} \quad x^2 + y^2 \leq r_l^2 l_0^2,\\
        0 & \mathrm{else},
    \end{cases}\\
\end{align*}$$

where $h_0 = 150 \ \mathrm{m}$, $l_0 = 5 \ \mathrm{km}$, $r_h = 2$, and $r_l = 10$. This is decomposed into a large-scale part $h_\mathrm{b}$ and a small-scale part with the spectral amplitude $h_\mathrm{w}$, such that

$$\begin{align*}
    h_\mathrm{b} \left(x, y\right) & = r_h h_\mathrm{w} \left(x, y\right),\\
    h_\mathrm{w} \left(x, y\right) & = \begin{cases}
        \frac{h_0}{2 \left(r_h + 1\right)} \left[1 + \cos \left(\frac{\pi}{r_l l_0} \sqrt{x^2 + y^2}\right)\right] & \mathrm{if} \quad x^2 + y^2 \leq r_l^2 l_0^2,\\
        0 & \mathrm{else}.
    \end{cases}
\end{align*}$$

The large-scale part is resolved, so that the grid is defined from it, whereas the small-scale part is used by MS-GWaM to parameterize the mountain waves generated by the resolved wind crossing it. As in the first mountain-wave example, the atmosphere is isothermal, with the default temperature $T_0 = 300 \ \mathrm{K}$ and the initial wind $\boldsymbol{u}_0 = \left(10, 0, 0\right)^\mathrm{T} \ \mathrm{m \ s^{- 1}}$.

The damping coefficient of the sponge is given by

$$\alpha_\mathrm{R} \left(x, y, z\right) = \frac{\alpha_{\mathrm{R}, \max}}{3} \left[\exp \left(\frac{\left|x\right| - L_x / 2}{\Delta x_\mathrm{R}}\right) + \exp \left(\frac{\left|y\right| - L_y / 2}{\Delta y_\mathrm{R}}\right) + \exp \left(\frac{z - L_z}{\Delta z_\mathrm{R}}\right)\right],$$

where $\alpha_{\mathrm{R}, \max} = 0.0179 \ \mathrm{s^{- 1}}$, $\Delta x_\mathrm{R} = L_x / 20$, $\Delta y_\mathrm{R} = L_y / 20$ and $\Delta z_\mathrm{R} = L_z / 10$. In contrast to the sinusoidal sponge discussed in the first example, this sponge applies a damping everywhere in the domain (weakest at the center of the surface, strongest in the upper corners). Once again, the sponge relaxes the wind to its initial state.

MS-GWaM is used with most of its parameters set to their default values. This means that the orographic source launches exactly one ray volume in each surface grid cell with a nonzero $h_\mathrm{w}$. Thus, the number of ray volumes allowed per grid cell (before merging is triggered) is `multiplication_factor` (a parameter of the WKB namelist) cubed, which is $4^3 = 64$.

Instead of a contour plot, the above script generates a scatter plot that visualizes the ray volumes, with the color representing the value of the phase-space wave-action density (see below).

![](examples/results/wkb_mountain_wave.svg)

## WKB wave packet

The script

```julia
# examples/scripts/wkb_wave_packet.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

x_size = 16
y_size = 16
z_size = 32

lx = 20000.0
ly = 20000.0
lz = 40000.0

rx = 0.25
ry = 0.25
rz = 0.25

x0 = 0.0
y0 = 0.0
z0 = 20000.0

a0 = 0.05

k = 16 * pi / lx
l = 16 * pi / ly
m = 32 * pi / lz

model = Compressible()
background = Realistic()
coriolis_frequency = 0.0001

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)
domain = DomainNamelist(;
    x_size,
    y_size,
    z_size,
    lx,
    ly,
    lz,
    base_comm = MPI.COMM_SELF,
)
auxiliary_state = State(Namelists(; atmosphere, domain))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

include("wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    save_ray_volumes = true,
    output_file = "wkb_wave_packet.h5",
    tmax = 900.0,
    output_interval = 900.0,
)
wkb = WKBNamelist(;
    wkb_mode = MultiColumn(),
    initial_wave_field = (alpha, x, y, z) ->
        (k, l, m, omega(x, y, z), wave_action_density(x, y, z)),
)

integrate(Namelists(; atmosphere, domain, output, wkb))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wkb_wave_packet.h5") do data
        plot_output(
            "examples/results/wkb_wave_packet.svg",
            data,
            ("nr", 8, 8, 16, 2);
            time_unit = "min",
        )
        return
    end
end

```

initializes an unresolved gravity-wave packet (i.e. one that is parameterized by MS-GWaM) in the stratosphere of a "realistic" compressible atmosphere (isentropic troposphere and isothermal stratosphere) and visualizes the resulting ray-volume distribution after fifteen minutes integration time (see below). Like the wave-packet script discussed above, it constructs an auxiliary state and uses helper functions to satisfy the gravity-wave dispersion and polarization relations.

![](examples/results/wkb_wave_packet.svg)

## Wave-packet helper functions

The script

```julia
# examples/scripts/wave_packet_tools.jl

function ijk(x, y, z)
    i = argmin(abs.(x .- auxiliary_state.grid.x .* lref))
    j = argmin(abs.(y .- auxiliary_state.grid.y .* lref))
    k = argmin(abs.(z .- auxiliary_state.grid.zc[i, j, :] .* lref))

    return CartesianIndex(i, j, k)
end

function rhobar(x, y, z)
    return auxiliary_state.atmosphere.rhobar[ijk(x, y, z)] .* rhoref
end

function thetabar(x, y, z)
    return auxiliary_state.atmosphere.thetabar[ijk(x, y, z)] .* thetaref
end

function n2(x, y, z)
    return auxiliary_state.atmosphere.n2[ijk(x, y, z)] ./ tref .^ 2
end

function envelope(x, y, z)
    r =
        sqrt(
            (rx * k * (x - x0))^2 +
            (ry * l * (y - y0))^2 +
            (rz * m * (z - z0))^2,
        ) / pi
    if r <= 1
        return (1 + cos(pi * r)) / 2
    else
        return 0.0
    end
end

function phi(x, y, z)
    return k * x + l * y + m * z
end

function omega(x, y, z)
    return -sqrt(
        (n2(x, y, z) * (k^2 + l^2) + coriolis_frequency^2 * m^2) /
        (k^2 + l^2 + m^2),
    )
end

function bhat(x, y, z)
    return a0 * n2(x, y, z) / m * envelope(x, y, z)
end

function uhat(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im / m / n2(x, y, z) * (omega(x, y, z)^2 - n2(x, y, z)) /
           (omega(x, y, z)^2 - coriolis_frequency^2) *
           (k * omega(x, y, z) + 1im * l * coriolis_frequency) *
           bhat(x, y, z)
end

function vhat(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im / m / n2(x, y, z) * (omega(x, y, z)^2 - n2(x, y, z)) /
           (omega(x, y, z)^2 - coriolis_frequency^2) *
           (l * omega(x, y, z) - 1im * k * coriolis_frequency) *
           bhat(x, y, z)
end

function what(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im * omega(x, y, z) / n2(x, y, z) * bhat(x, y, z)
end

function pihat(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           kappa / rsp / thetabar(x, y, z) * 1im / m *
           (omega(x, y, z)^2 - n2(x, y, z)) / n2(x, y, z) * bhat(x, y, z)
end

function wave_action_density(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           rhobar(x, y, z) / 2 * omega(x, y, z) * (k^2 + l^2 + m^2) /
           n2(x, y, z)^2 / (k^2 + l^2) * bhat(x, y, z)^2
end

```

provides helper functions that implement the gravity-wave dispersion and polarization relations needed for the initialization of wave packets. It is to be included below the construction of a corresponding auxiliary state and the extraction of $\kappa = R / c_p$, $R$ and $g$ from its `Constants` instance.
