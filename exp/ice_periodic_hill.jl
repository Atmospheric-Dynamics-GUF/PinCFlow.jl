# examples/scripts/periodic_hill.jl

using Pkg

#Pkg.activate("examples")

using MPI
using HDF5
#using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

h0 = 500.0
l0 = 10000.0

lz = 20000.0
zr = 10000.0

atmosphere = AtmosphereNamelist(;
    model = Boussinesq(),
    background = StableStratification(),
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
)
domain = DomainNamelist(; x_size = 40, z_size = 40, lx = 20000.0, lz, npx, npz)
grid = GridNamelist(;
    resolved_topography = (x, y) -> h0 / 2 * (1 + cos(pi / l0 * x)),
)
output =
    OutputNamelist(; output_variables = (:w,), output_file = "periodic_hill.h5")
sponge = SpongeNamelist(;
    rhs_sponge = (x, y, z, t, dt) ->
        z >= zr ? sin(pi / 2 * (z - zr) / (lz - zr))^2 / dt : 0.0,
)
discretization = DiscretizationNamelist(;
    dtmax = 10.0,
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge, discretization))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("periodic_hill.h5") do data
        plot_output(
            "examples/results/periodic_hill.svg",
            data,
            ("w", 1, 1, 1, 2);
        )
        return
    end
end
