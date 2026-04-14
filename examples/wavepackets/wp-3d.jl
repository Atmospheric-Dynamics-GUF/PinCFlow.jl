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

x_size = 512
y_size = 16
z_size = 1600

lx = 9000e3
ly = 300e3
lz = 100e3

rx = 1500e3
ry = 0
rz = 5e3

x0 = 0.0
y0 = 0.0
z0 = 30e3

a0 = 1.0

k = 0
l = 2 * pi / 300e3
m = 2 * pi / 3e3

dzr = lz / 10
alpharmax = 0.0179

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

atmosphere = AtmosphereNamelist(;
    background,
    model,
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
    output_variables = [:u, :v, :w, :rhop],
    output_file = "wp-3d.h5",
    tmax = 1000.0,
    output_interval = 1000.0,
)

discretization = DiscretizationNamelist(; dtmax = 100)

turbulence = TurbulenceNamelist(;
    turbulence_scheme = :NoTurbulence,
    initial_tke = (x, y, z) -> qtilde(x, y, z) / 2,
)

sponge = SpongeNamelist(;
    lhs_sponge = (x, y, z, t, dt) -> alpharmax * exp((z - lz) / dzr),
)

tracer = TracerNamelist(;
    tracer_setup = :TracerOn,
    leading_order_impact = true,
    initial_tracer = (x, y, z) ->
        real(chihat(x, y, z) * exp(1im * phi(x, y, z))) + z,
    background_tracer = (x, y, z) -> z,
)

integrate(
    Namelists(;
        atmosphere,
        domain,
        output,
        tracer,
        sponge,
        turbulence,
        discretization,
    ),
)
