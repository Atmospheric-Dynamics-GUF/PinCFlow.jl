# examples/submit/mountain_wave.jl

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
    ndy = 40,
    ndz = 40,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
    npx = 8,
    npy = 8,
)
grid = GridNamelist(; mountain_case = 4)
output = OutputNamelist(; output_variables = (:w,), output_file = output_file)
sponge = SpongeNamelist(;
    spongelayer = true,
    alpharmax = 1.79E-2,
    betarmax = 0.0E+0,
    lateralsponge = true,
    spongetype = SinusoidalSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
