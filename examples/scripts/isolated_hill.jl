# examples/scripts/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

h0 = 500.0
l0 = 10000.0

lx = 80000.0
lz = 30000.0
zr = 20000.0

atmosphere = AtmosphereNamelist(;
    model = PseudoIncompressible(),
    background = LapseRates(),
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
)
domain = DomainNamelist(; x_size = 40, z_size = 120, lx, lz, npx, npz)
grid = GridNamelist(; resolved_topography = (x, y) -> h0 / (1 + x^2 / l0^2))
output =
    OutputNamelist(; output_variables = (:w,), output_file = "isolated_hill.h5")
sponge = SpongeNamelist(;
    rhs_sponge = (x, y, z, t, dt) ->
        z >= zr ? sin(pi / 2 * (z - zr) / (lz - zr))^2 / dt : 0.0,
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("isolated_hill.h5") do data
        plot_output(
            "examples/results/isolated_hill.svg",
            data,
            ("w", 1, 1, 1, 2);
        )
        return
    end
end
