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

tau_q_sink = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 0.0
@info "tau_q_sink" tau_q_sink

period_ss = length(ARGS) ≥ 6 ? parse(Float64, ARGS[6]) : 3600.0
@info "period_ss" period_ss

run = length(ARGS) ≥ 7 ? ARGS[7] : "0000_00"
@info "run" run

#outfile = "/home/b/b383844/PinCFlow/sedimentation/results/ice_mountain_wave_$(run).h5"
outfile = "/work/bb1097/b383844/PinCFlow/qv_forcing/results/ice_mountain_wave_$(run).h5"

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
	icesetup = IceOn(),
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
	tmax = 2000.0,
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
    relaxed_u = (x, y, z, t, dt) -> 10.0 * cos(2 * pi * t/ period_ss),
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
