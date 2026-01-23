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

x_size = 32 #x_size = 16
y_size = 1 #y_size = 16
z_size = 64 #z_size = 32


lx = 20000.0
ly = 20000.0
lz = 40000.0

rx = 0.25
ry = 0.25
rz = 0.25

x0 = 0.0
y0 = 0.0
z0 = 20000.0

x_c = 0.0
y_c = 0.0
z_c = 20000.0

sigma_x = lx / 6
sigma_y = ly / 6
sigma_z = lz / 12

a0 = 0.05

k = 16 * pi / lx
l = 16 * pi / ly
m = 32 * pi / lz

model = Boussinesq()
background = StableStratification()
coriolis_frequency = 0.0

atmosphere = AtmosphereNamelist(; model, background, coriolis_frequency)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
auxiliary_state = State(Namelists(; atmosphere, domain))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

#println(auxiliary_state.grid.kp)
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
output = OutputNamelist(;
    output_variables = (:u, :v, :w),
    output_file = "wave_packet.h5",
    save_ray_volumes = true,
    tmax = 100.0,
    output_interval = 20.0,
)

integrate(Namelists(; atmosphere, domain, output))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wave_packet.h5") do data
        plot_output(
            "examples/results/wave_packet.svg",
            data,
            ("w", 4, 1, 8, 1),
            #("v", 20, 20, 40, 2),
            #("w", 20, 20, 40, 2);
            time_unit = "min",
        )
        return
    end
end
