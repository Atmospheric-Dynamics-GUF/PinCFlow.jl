# examples/scripts/wp-3d.jl

using Pkg

using MPI
using HDF5
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

x_size = 1
y_size = 1
z_size = 30

lx = 1e3
ly = 1e3
lz = 30e3

rx = 0.0
ry = 0.0
rz = 2e3

x0 = 0.0
y0 = 0.0
z0 = 10e3

a0 = 0.9

k = 2 * pi / 1e3
l = 0
m = 2 * pi / 1e3

background = :Isothermal
model = :Compressible
coriolis_frequency = 0

atmosphere = AtmosphereNamelist(; model, coriolis_frequency, background)
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
    output_variables = [:u, :v, :w, :rhop, :e, :dtkedt, :dchidt],
    output_file = "wkb-wp-1d.h5",
    tmax = 3600,
    output_interval = 60,
)

wkb = WKBNamelist(;
    use_saturation = false,
    nrz = 80,
    wkb_mode = :MultiColumn,
    filter_order = 2,
    filter_type = :BoxFilter,
    initial_wave_field = (alpha, x, y, z) ->
        (k, l, m, omega(x, y, z), wave_action_density(x, y, z)),
    turbulence_damping = false,
)
turbulence = TurbulenceNamelist(;
    turbulence_scheme = :NoTurbulence,
    tracer_coupling = true,
)

tracer = TracerNamelist(;
    tracer_setup = :TracerOn,
    leading_order_impact = false,
    next_order_impact = true,
    turbulence_impact = false,
    initial_tracer = (x, y, z) -> z,
)

integrate(Namelists(;
    atmosphere,
    domain,
    output,
    wkb,
    tracer,
    turbulence,
))