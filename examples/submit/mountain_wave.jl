# examples/submit/mountain_wave.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

atmosphere = AtmosphereNamelist(; initial_wind = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 40,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
    npx = 8,
    npy = 8,
)
grid = GridNamelist(; mountain_case = 4)
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
