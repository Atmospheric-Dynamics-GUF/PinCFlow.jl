# examples/scripts/ice_mountain_2D.jl
# Run with:
# julia --project examples/scripts/ice_mountain_2D.jl

using Pkg

#Pkg.activate("examples")

using MPI
using HDF5
#using CairoMakie
#using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

h0 = 1000.0
l0 = 1000.0

lx = 40000.0
ly = 20000.0
lz = 20000.0
dxr = lx / 2
dyr = ly / 2
dzr = lz / 2
alpharmax = 0.0179

MPI.Init()
rank = MPI.Comm_rank(MPI.COMM_WORLD)
base_comm = MPI.Comm_split(MPI.COMM_WORLD, rank % 2, rank)

if rank % 2 == 0
    sigma = 1
    output_file = "ice_mountain_wave_positive_u.h5"
else
    sigma = -1
    output_file = "ice_mountain_wave_negative_u.h5"
end

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0 * sigma,
)
domain = DomainNamelist(;
    x_size = 40,
    y_size = 1,
    z_size = 80,
    lx,
    ly,
    lz,
    npx,
    npy,
    npz,
    base_comm
)
grid = GridNamelist(;
    resolved_topography = (x, y) -> h0 / (1 + x^2 / l0^2),
)
ice = IceNamelist(;
	icesetup = IceOn(),
#	ice_test_case = MultipleWavePackets(),
	dt_ice = 2.0,
	nscx = 1,
	nscy = 1,
	nscz = 1,
	cloudcover = CloudCoverOff(),
)
output =
    OutputNamelist(; output_variables = (:w, :u, :n, :qv, :q, :iaux1, :iaux2, :iaux3), 
    output_steps = false,
	output_interval = 100.0,
	tmax = 2000.0,
    output_file)

    sponge = SpongeNamelist(;
    lhs_sponge = (x, y, z, t, dt) -> begin
        alpharx =
            abs(x) >= (lx - dxr) / 2 ?
            sin(pi * (abs(x) - (lx - dxr) / 2) / dxr)^2 : 0.0
        alphary = 0.0
#            abs(y) >= (ly - dyr) / 2 ?
#            sin(pi * (abs(y) - (ly - dyr) / 2) / dyr)^2 : 0.0
        alpharz =
            z >= lz - dzr ? sin(pi / 2 * (z - (lz - dzr)) / dzr)^2 : 0.0
        return alpharmax * (alpharx + alphary + alpharz) / 3
    end,
    relaxed_u = (x, y, z, t, dt) -> sigma * 10.0 * cos(pi * t/ 1800.0),
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge, ice))

MPI.Barrier(MPI.COMM_WORLD)

# if MPI.Comm_rank(MPI.COMM_WORLD) == 0
#     h5open("mountain_wave.h5") do data
#         plot_output(
#             "examples/results/mountain_wave.svg",
#             data,
#             ("w", 20, 20, 10, 2);
#         )
#         return
#     end
# end
