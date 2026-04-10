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

x_size = 32
y_size = 32
z_size = 1066

lx = 30e3
ly = 30e3
lz = 100e3

rx = 0.0
ry = 0.0
rz = 0.05

x0 = 0.0
y0 = 0.0
z0 = 20e3

a0 = 0.9

k = 2 * pi / 30e3
l = 2 * pi / 30e3
m = 2 * pi / 3e3

background = :Isothermal
model = :Compressible
coriolis_frequency = 1e-4

atmosphere = AtmosphereNamelist(; background, coriolis_frequency)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
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

turbulence = TurbulenceNamelist(; turbulence_scheme = TKEScheme())

output = OutputNamelist(;
    output_variables = (:u, :v, :w, :rhop),
    output_file = "wp-1d-turbulence.h5",
    tmax = 3600 * 5,
    output_interval = 360,
)

tracer = TracerNamelist(;
    tracer_setup = TracerOn(),
    initial_tracer = (x, y, z) ->
        real(chihat(x, y, z) * exp(1im * phi(x, y, z))) + z,
)

integrate(
    Namelists(; atmosphere, domain, output, tracer, turbulence),
)
