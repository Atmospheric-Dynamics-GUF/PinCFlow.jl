# examples/submit/wkb_mountain_wave.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

atmosphere = AtmosphereNamelist(; initial_wind = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 40,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx = 4.0E+5,
    ly = 4.0E+5,
    lz = 2.0E+4,
    npx = 8,
    npy = 8,
)
grid = GridNamelist(;
    mountain_height = 1.5E+2,
    mountain_half_width = 5.0E+3,
    mountain_case = 13,
    height_factor = 2.0E+0,
    width_factor = 1.0E+1,
)
output = OutputNamelist(;
    output_variables = (:w,),
    output_file = "wkb_mountain_wave.h5",
)
setting = SettingNamelist(; test_case = WKBMountainWave())
sponge = SpongeNamelist(;
    use_sponge = true,
    sponge_extent = 1.0E-1,
    alpharmax = 1.79E-2,
    betarmax = 0.0E+0,
    lateral_sponge = true,
    sponge_type = ExponentialSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)
wkb = WKBNamelist(;
    xrmin = -2.0E+5,
    xrmax = 2.0E+5,
    yrmin = -2.0E+5,
    yrmax = 2.0E+5,
    zrmin = 0.0,
    zrmax = 2.0E+4,
)

integrate(Namelists(; atmosphere, domain, grid, output, setting, sponge, wkb))
