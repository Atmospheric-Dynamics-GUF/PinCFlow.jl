# examples/scripts/wave_packet_tke.jl

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

x_size = 128
y_size = 16
z_size = 160

lx = 900000.0
ly = 100000.0
lz = 30000.0

rx = 0.1
ry = 0.0
rz = 0.1

x0 = 0.0
y0 = 0.0
z0 = 10e3

a0 = 1

k = 2 * pi / 100e3
l = 2 * pi / 100e3
m = 2 * pi / 3e3
background = Isothermal()
model = PseudoIncompressible()
coriolis_frequency = 1.e-4

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

include("wave_packet_tools.jl")

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
    output_variables = (:u, :v, :w, :chi),
    output_file = "wp-test.h5",
    tmax = 36000.0,
    output_interval = 36000.0,
)
tracer = TracerNamelist(;
    tracer_setup = TracerOn(),
    initial_tracer = (x, y, z) ->
        real(chihat(x, y, z) * exp(1im * phi(x, y, z))) + z,
)
poisson = PoissonNamelist(; tolerance = 1.e-6, poisson_iterations = 1000)
turbulence = TurbulenceNamelist(;
    turbulence_scheme = TKEScheme(),
    momentum_coupling = true,
)

integrate(Namelists(; atmosphere, domain, output, tracer, turbulence, poisson))