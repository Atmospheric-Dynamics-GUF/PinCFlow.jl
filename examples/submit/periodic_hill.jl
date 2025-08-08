# examples/submit/periodic_hill.jl

using PinCFlow

atmosphere = AtmosphereNamelist(; backgroundflow_dim = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
    sizex = 40,
    sizey = 1,
    sizez = 40,
    lx_dim = (-1.0E+4, 1.0E+4),
    ly_dim = (-1.0E+4, 1.0E+4),
    lz_dim = (0.0E+0, 2.0E+4),
)
grid = GridNamelist(; mountainheight_dim = 1.0E+1, mountainwidth_dim = 1.0E+4)
output = OutputNamelist(;
    output_variables = (:w,),
    output_file = ARGS[1] * "/pincflow_output.h5",
)
sponge = SpongeNamelist(; spongelayer = true)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
