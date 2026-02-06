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

x_size = 1
y_size = 1
z_size = 20

lx = 30e3
ly = 30e3
lz = 60e3

rx = 0.0
ry = 0.0
rz = 0.05

x0 = 0.0
y0 = 0.0
z0 = 20e3

a0 = 2

k = 2 * pi / 30e3
l = 2 * pi / 30e3
m = 2 * pi / 3e3

background = Isothermal()
model = PseudoIncompressible()
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

include("wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

output = OutputNamelist(;
    save_ray_volumes = true,
    output_variables = (:u, :v, :w, :rhop, :dtkedt),
    output_file = "wkb-wp-3d-nobp.h5",
    tmax = 360.0,
    output_interval = 36.0,
)

discretization = DiscretizationNamelist(; dtmax = 1.0)

wkb = WKBNamelist(;
    nrz = 1,
    use_saturation = false,
    wkb_mode = MultiColumn(),
    initial_wave_field = (alpha, x, y, z) ->
        (k, l, m, omega(x, y, z), wave_action_density(x, y, z)),
    filter_order = 3,
)

turbulence = TurbulenceNamelist(;
    turbulence_scheme = TKEScheme(),
    momentum_coupling = true,
    wave_action_coupling = true,
    initial_tke = (x, y, z) -> 5e-5,
)

integrate(
    Namelists(; atmosphere, domain, output, turbulence, discretization, wkb),
)
