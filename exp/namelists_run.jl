# examples/scripts/mountain_wave.jl

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

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
)
domain = DomainNamelist(;
    x_size = 4,
    y_size = 1,
    z_size = 8,
    lx,
    ly,
    lz,
    npx,
    npy,
    npz,
)
grid = GridNamelist(;
    resolved_topography = (x, y) -> h0 / (1 + x^2 / l0^2),
)
ice = IceNamelist(;
	ice_setup = IceOn(),
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
    output_file = "ice_mountain_wave.h5")

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
    relaxed_u = (x, y, z, t, dt) -> 10.0 * cos(pi * t/ 1800.0),
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge, ice))

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

tau_q_sink = 0.00000000001

tau_qv_source = 3000.0

run = 2901_06

tmax = 2.0e5

#outfile = "/home/b/b383844/PinCFlow/sedimentation/results/ice_mountain_wave_$(run).h5"
outfile = "/work/bb1097/b383844/PinCFlow/adv/results/ice_mountain_wave_$(run).h5"

h0 = 1000.0
l0 = 1000.0

lx = 40000.0
ly = 20000.0
lz = 20000.0
dxr = lx / 2
dyr = ly / 2
dzr = lz / 2
alpharmax = 0.0179

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
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
)
grid = GridNamelist(;
    resolved_topography = (x, y) -> h0 / (1 + x^2 / l0^2),
)
ice = IceNamelist(;
    tau_q_sink = tau_q_sink,
    tau_qv_source = tau_qv_source,
	ice_setup = IceOn(),
#	ice_test_case = MultipleWavePackets(),
	dt_ice = 2.0,
	nscx = 1,
	nscy = 1,
	nscz = 1,
	cloudcover = CloudCoverOff(),
)
output =
    OutputNamelist(; output_variables = (:w, :u, :n, :qv, :q, :iaux1, :iaux2, :iaux3, :clc), 
    output_steps = false,
	output_interval = 100.0,
	tmax = tmax,
    output_file = outfile)

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
    relaxed_u = (x, y, z, t, dt) -> 10.0 + 10.0 * sin(2 * pi * t/ 1.0e5),
)

# save sbatch script copy and ice_mountain_2D.jl to output directory
MPI.Init()
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    output_dir = dirname(outfile)
    sbatch_src = "/home/b/b383844/PinCFlow/PinCFlow.jl/examples/scripts/levante/ice_mountain_2D_sinks.sh"
    sbatch_dst = "/work/bb1097/b383844/PinCFlow/batch/ice_mountain_2D_$(run).sbatch"
    cp(sbatch_src, sbatch_dst; force=true)
    script_src = "/home/b/b383844/PinCFlow/PinCFlow.jl/exp/ice_mountain_2D.jl"
    script_dst = "/work/bb1097/b383844/PinCFlow/julia/ice_mountain_2D_$(run).jl"
    cp(script_src, script_dst; force=true)
end

integrate(Namelists(; atmosphere, domain, grid, output, sponge, ice))

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

