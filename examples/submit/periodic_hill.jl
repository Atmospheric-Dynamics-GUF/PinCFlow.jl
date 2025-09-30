# examples/submit/periodic_hill.jl

using PinCFlow

@ivy if length(ARGS) == 0
    output_file = "./pincflow_output.h5"
elseif length(ARGS) == 1
    output_file = ARGS[1] * "/pincflow_output.h5"
else
    error("Too many arguments to the script!")
end

atmosphere = AtmosphereNamelist(; initial_wind = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
    ndx = 40,
    ndy = 1,
    ndz = 40,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
)
grid = GridNamelist(; mountainheight_dim = 1.0E+1, mountainwidth_dim = 1.0E+4)
output = OutputNamelist(; output_variables = (:w,), output_file = output_file)
sponge = SpongeNamelist(; spongelayer = true)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
