# examples/submit/mountain_wave.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

atmosphere = AtmosphereNamelist(;
    initial_wind = (1.0E+1, 0.0E+0, 0.0E+0),
    coriolis_frequency = 0.0E+0,
)
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 40,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
    npx,
    npy,
    npz,
)
grid = GridNamelist(;
    resolved_topography = (namelists, x, y) ->
        100 / (1 + (x^2 + y^2) / 1000^2),
)
output =
    OutputNamelist(; output_variables = (:w,), output_file = "mountain_wave.h5")
sponge = SpongeNamelist(;
    alpharmax = 1.79E-2,
    lateral_sponge = true,
    sponge_type = SinusoidalSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
