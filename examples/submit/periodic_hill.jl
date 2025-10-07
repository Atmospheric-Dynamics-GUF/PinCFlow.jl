# examples/submit/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

atmosphere = AtmosphereNamelist(;
    initial_wind = (1.0E+1, 0.0E+0, 0.0E+0),
    coriolis_frequency = 0.0E+0,
)
domain = DomainNamelist(;
    x_size = 40,
    y_size = 1,
    z_size = 40,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
    npx,
    npz,
)
grid = GridNamelist(; mountain_height = 1.0E+1, mountain_half_width = 1.0E+4)
output =
    OutputNamelist(; output_variables = (:w,), output_file = "periodic_hill.h5")
sponge = SpongeNamelist(; betarmax = 1.0E+0)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
