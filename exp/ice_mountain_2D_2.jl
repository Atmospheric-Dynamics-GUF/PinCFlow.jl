# examples/scripts/wkb_mountain_wave.jl

using Pkg

using MPI
using HDF5
#using CairoMakie
#using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

h0 = 150.0
l0 = 5000.0
rl = 10
rh = 2

lx = 4.0e5
ly = 20000.0
lz = 20000.0
dxr = lx / 20
dyr = ly / 20
dzr = lz / 10
alpharmax = 0.0179

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
)
domain = DomainNamelist(;
    x_size = 40,
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
        h0 / 2 * (1 + cos(pi / (rl * l0) * abs(x))) : 0.0,
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
output = OutputNamelist(; 
    output_variables = (:w, :u, :n, :qv, :q, :iaux1, :iaux2, :iaux3), 
    output_steps = false,
	output_interval = 100.0,
	tmax = 2000.0,
    save_ray_volumes = true,
    output_file = "ice_mountain_wave.h5",
)
sponge = SpongeNamelist(;
    lhs_sponge = (x, y, z, t, dt) ->
        alpharmax / 3 * (
            exp((abs(x) - lx / 2) / dxr) +
         #   exp((abs(y) - ly / 2) / dyr) +
            exp((z - lz) / dzr)
        ),
    relaxed_u = (x, y, z, t, dt) -> 10.0,
)
#wkb = WKBNamelist(; wkb_mode = MultiColumn())

integrate(Namelists(; atmosphere, domain, grid, output, sponge))

# if MPI.Comm_rank(MPI.COMM_WORLD) == 0
#     h5open("wkb_mountain_wave.h5") do data
#         plot_output(
#             "examples/results/wkb_mountain_wave.svg",
#             data,
#             ("nr", 20, 20, 10, 2);
#         )
#         return
#     end
# end
