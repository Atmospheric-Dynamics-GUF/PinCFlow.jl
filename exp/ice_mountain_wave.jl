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

h0 = 100.0
l0 = 1000.0

lx = 20000.0
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
    y_size = 40,
    z_size = 40,
    lx,
    ly,
    lz,
    npx,
    npy,
    npz,
)
grid = GridNamelist(;
    resolved_topography = (x, y) -> h0 / (1 + (x^2 + y^2) / l0^2),
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
    OutputNamelist(; output_variables = (:w, :n, :qv, :q, :iaux1, :iaux2, :iaux3), 
    output_steps = false,
	output_interval = 50.0,
	tmax = 100.0,
    output_file = "exp/mountain_wave.h5")

    sponge = SpongeNamelist(;
    lhs_sponge = (x, y, z, t, dt) -> begin
        alpharx =
            abs(x) >= (lx - dxr) / 2 ?
            sin(pi * (abs(x) - (lx - dxr) / 2) / dxr)^2 : 0.0
        alphary =
            abs(y) >= (ly - dyr) / 2 ?
            sin(pi * (abs(y) - (ly - dyr) / 2) / dyr)^2 : 0.0
        alpharz =
            z >= lz - dzr ? sin(pi / 2 * (z - (lz - dzr)) / dzr)^2 : 0.0
        return alpharmax * (alpharx + alphary + alpharz) / 3
    end,
    relaxed_u = (x, y, z, t, dt) -> 10.0,
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge, ice))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("mountain_wave.h5") do data
        plot_output(
            "examples/results/mountain_wave.svg",
            data,
            ("w", 20, 20, 10, 2);
        )
        return
    end
end
