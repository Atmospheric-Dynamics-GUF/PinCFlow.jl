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

x_size = 16 #x_size = 16
y_size = 1 #y_size = 16
z_size = 32 #z_size = 32

lx = 20000.0
ly = 20000.0
lz = 40000.0

rx = 0.25
ry = 0.25
rz = 0.25

x0 = 0.0
y0 = 0.0
z0 = 10000.0

a0 = 0.05

k = 16 * pi / lx
l = 16 * pi / ly
m = 32 * pi / lz

model = Boussinesq()
background = StableStratification()
coriolis_frequency = 0.0001

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
auxiliary_state = State(Namelists(; atmosphere, domain))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

include("wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)
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
            ("nr", 4, 1, 8, 1);
            time_unit = "min",
        )
        return
    end
end
