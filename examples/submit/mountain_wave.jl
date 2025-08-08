# examples/submit/mountain_wave.jl

using PinCFlow

atmosphere = AtmosphereNamelist(; backgroundflow_dim = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
    sizex = 40,
    sizey = 40,
    sizez = 40,
    lx_dim = (-1.0E+4, 1.0E+4),
    ly_dim = (-1.0E+4, 1.0E+4),
    lz_dim = (0.0, 2.0E+4),
    npx = 8,
    npy = 8,
)
grid = GridNamelist(; mountain_case = 4)
output = OutputNamelist(;
    output_variables = (:w,),
    output_file = ARGS[1] * "/pincflow_output.h5",
)
sponge = SpongeNamelist(;
    spongelayer = true,
    spongealphaz_dim = 1.79E-2,
    unifiedsponge = true,
    lateralsponge = true,
    spongetype = SinusoidalSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
