# examples/scripts/wp-3d.jl

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

x_size = 52
y_size = 1
z_size = 100

lx = 9000e3
ly = 300e3
lz = 100e3

rx = 1500e3
ry = 0
rz = 5e3

x0 = 0.0
y0 = 0.0
z0 = 30e3

a0 = 0.5

k = 0
l = 2 * pi / 300e3
m = 2 * pi / 3e3

background = :Isothermal
model = :Compressible
coriolis_frequency = 1e-4

atmosphere = AtmosphereNamelist(; background, coriolis_frequency, model)
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
(; lturb) = auxiliary_state.turbulence.turbulenceconstants

include("wave_packet_tools-3d.jl")

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

output = OutputNamelist(;
    save_ray_volumes = false,
    output_variables = (:u, :v, :w, :rhop, :dchidt, :e),
    output_file = "wkb-wp-3d.h5",
    tmax = 1000.0,
    output_interval = 1000.0,
)

wkb = WKBNamelist(;
    use_saturation = false,
    nrz = 1,
    wkb_mode = MultiColumn(),
    initial_wave_field = (alpha, x, y, z) ->
        (k, l, m, omega(x, y, z), wave_action_density(x, y, z)),
    turbulence_damping = true,
)

discretization = DiscretizationNamelist(; dtmax = 100)

turbulence = TurbulenceNamelist(; turbulence_scheme = TKEScheme())

tracer = TracerNamelist(;
    tracer_setup = TracerOn(),
    leading_order_impact = true,
    next_order_impact = false,
    turbulence_impact = false,
    initial_tracer = (x, y, z) -> z,
)

integrate(
    Namelists(;
        atmosphere,
        domain,
        output,
        wkb,
        tracer,
        turbulence,
        discretization,
    ),
)
