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

run = "2705_01"

#outfile = "/home/b/b383844/PinCFlow/sedimentation/results/mountain_wave_$(run).h5"
outfile = "/work/bb1097/b383844/PinCFlow/adv/results/mountain_wave_$(run).h5"

tmax = 2.0e3

h0 = 150.0
l0 = 5000.0
rl = 10
rh = 2

lx = 400_000.0
ly = 400_000.0
lz = 20_000.0
dxr = lx / 20
dyr = ly / 20
dzr = lz / 10
alpharmax = 0.0179



atmosphere = AtmosphereNamelist(;
    background = LapseRates(),
    temperature = 280.0,
    potential_temperature = 280.0,
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
)

domain = DomainNamelist(;
    x_size = 400,
    y_size = 1,
    z_size = 40,
    lx,
    ly,
    lz,
    npx,
    npy,
    npz,
)
grid = GridNamelist(;
    resolved_topography = (x, y) ->
        x^2 <= (rl * l0)^2 ?
        h0 / 2 * (1 + cos(pi / (rl * l0) * abs(x) ) ) * rh / (rh + 1) + 
        h0 / 2 * (1 + cos(pi / (rl * l0) * abs(x) ) ) / 
        (rh + 1) * cos(pi / l0 * abs(x) ) : 0.0,
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
output = OutputNamelist(; 
    output_variables = (:w, :u, :n, :nNuc, :qv, :q, :thetap, :pip, :iaux1, :iaux2, :iaux3, :iaux4, :iaux5, :clc), 
    output_steps = false,
	output_interval = 10.0,
	tmax = tmax,
    save_ray_volumes = true,
    output_file = outfile,
)

sponge = SpongeNamelist(;
    lhs_sponge = (x, y, z, t, dt) ->
        alpharmax / 3 * (
            exp((abs(x) - lx / 2) / dxr) +
 #           exp((abs(y) - ly / 2) / dyr) +
            exp((z - lz) / dzr)
        ),
    relaxed_u = (x, y, z, t, dt) -> 10.0 + (10.0 * sin(2 * pi * t/ 1.0e5)),
)

# save sbatch script copy and ice_mountain_2D.jl to output directory
MPI.Init()
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    output_dir = dirname(outfile)
    sbatch_src = "/home/b/b383844/PinCFlow/PinCFlow.jl/examples/scripts/levante/ice_dump.sh"
    sbatch_dst = "/work/bb1097/b383844/PinCFlow/batch/ice_dump_$(run).sbatch"
    cp(sbatch_src, sbatch_dst; force=true)
    script_src = "/home/b/b383844/PinCFlow/PinCFlow.jl/examples/scripts/mountain_wave.jl"
    script_dst = "/work/bb1097/b383844/PinCFlow/julia/mountain_wave_$(run).jl"
    cp(script_src, script_dst; force=true)
end

integrate(Namelists(; atmosphere, domain, grid, output, sponge, ice))

