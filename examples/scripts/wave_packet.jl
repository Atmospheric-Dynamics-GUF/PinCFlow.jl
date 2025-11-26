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
