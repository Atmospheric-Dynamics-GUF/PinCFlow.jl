# examples/submit/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0E+0,
    initial_u = (x, y, z) -> 1.0E+1,
)
domain = DomainNamelist(;
    x_size = 20,
    y_size = 1,
    z_size = 20,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
    npx,
    npz,
)
grid = GridNamelist(;
    resolved_topography = (x, y) -> 5 * (1 + cos(pi / 10000 * x)),
)
output =
    OutputNamelist(; output_variables = (:w,), output_file = "periodic_hill.h5")
sponge = SpongeNamelist(; betarmax = 1.0E+0)

integrate(Namelists(; atmosphere, domain, grid, output, sponge, turbulence, setting))
