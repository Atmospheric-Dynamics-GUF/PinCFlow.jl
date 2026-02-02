# examples/scripts/wave_packet_STIH.jl

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
y_size = 1
z_size = 640

lx = 100e3
ly = 100e3
lz = 100000.0

rx = 0.0
ry = 0.0
rz = 0.05

x0 = 0.0
y0 = 0.0
z0 = 30e3

a0 = 2

k = 2 * pi / 100e3
l = 0.0
m = 2 * pi / 5e3
background = Isothermal()
model = PseudoIncompressible()
coriolis_frequency = 0.0

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
    output_variables = (:u, :v, :w, :rhop),
    output_file = "wp-5km-nocoup.h5",
    tmax = 3600.0,
    output_interval = 36.0,
)
turbulence = TurbulenceNamelist(;
    turbulence_scheme = TKEScheme(),
    momentum_coupling = false,
    initial_tke = (x, y, z) -> qtilde(x, y, z) / 2,
)

integrate(Namelists(; atmosphere, domain, output, turbulence))