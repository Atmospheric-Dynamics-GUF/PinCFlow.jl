# src/Examples/wkb_mountain_wave.jl
# examples/scripts/wp-3d.jl

using Pkg

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
z_size = 40

h0 = 1000.0
l0 = 50000.0
rl = 10
rh = 2

lx = 4000000.0
ly = 4000000.0
lz = 20000.0
dxr = lx / 20
dyr = ly / 20
dzr = lz / 10
alpharmax = 0.0179

output_file = "wkb-mw-3d.h5"
visualize = false

atmosphere = AtmosphereNamelist(;
    background = :Isothermal,
    coriolis_frequency = 1e-4,
    initial_u = (x, y, z) -> 15.0,
    initial_v = (x, y, z) -> 10.0,
)

domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

grid = GridNamelist(;
    resolved_topography = (x, y) ->
        x^2 + y^2 <= (rl * l0)^2 ?
        h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) * rh / (rh + 1) : 0.0,
    unresolved_topography = (alpha, x, y) ->
        x^2 + y^2 <= (rl * l0)^2 ?
        (
            pi / l0,
            pi / l0,
            h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) / (rh + 1) * cos(pi * x / l0),
        ) : (0.0, 0.0, 0.0),
)

output = OutputNamelist(;
    output_file,
    output_variables = [:u, :v, :w, :dchidt, :e, :dtkedt],
    save_ray_volumes = false,
    tmax = 3.6e4 * 1.5,
    output_interval = 3.6e3,
)

sponge = SpongeNamelist(;
    lhs_sponge = (x, y, z, t, dt) ->
        alpharmax / 3 * (
            exp((abs(x) - lx / 2) / dxr) +
            exp((abs(y) - ly / 2) / dyr) +
            exp((z - lz) / dzr)
        ),
    relaxed_u = (x, y, z, t, dt) -> 10.0,
    relaxed_v = (x, y, z, t, dt) -> 10.0,
)

tracer = TracerNamelist(;
    tracer_setup = :TracerOn,
    next_order_impact = false,
    turbulence_impact = true,
    apply_sponge_to_tracer = true,
    initial_tracer = (x, y, z) -> z,
    background_tracer = (x, y, z) -> z,
)

wkb = WKBNamelist(;
    wkb_mode = :MultiColumn,
    turbulence_damping = true,
    use_saturation = false,
)

turbulence = TurbulenceNamelist(; turbulence_scheme = :TKEScheme)

integrate(
    Namelists(;
        atmosphere,
        domain,
        grid,
        output,
        sponge,
        wkb,
        tracer,
        turbulence,
    ),
)

if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open(output_file) do data
        plot_output("wkb-mv00.svg", data, ("nr", 40, 20, 15, 1))
        plot_output("wkb-mv02.svg", data, ("nr", 40, 20, 15, 3))
        plot_output("wkb-mv04.svg", data, ("nr", 40, 20, 15, 5))
        plot_output("wkb-mv06.svg", data, ("nr", 40, 20, 15, 7))
        plot_output("wkb-mv08.svg", data, ("nr", 40, 20, 15, 9))
        plot_output("wkb-mv10.svg", data, ("nr", 40, 20, 15, 11))
        return
    end
end