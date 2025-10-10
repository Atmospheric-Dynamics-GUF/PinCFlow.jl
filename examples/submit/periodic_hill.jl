# examples/submit/periodic_hill.jl

using PinCFlow

@ivy if length(ARGS) == 0
    output_file = "./pincflow_output.h5"
elseif length(ARGS) == 1
    output_file = ARGS[1] * "/pincflow_output.h5"
else
    error("Too many arguments to the script!")
end

atmosphere = AtmosphereNamelist(; background = UniformBoussinesq())
domain = DomainNamelist(;
    x_size = 20,
    y_size = 1,
    z_size = 20,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
)
setting = SettingNamelist(; model = Boussinesq())
grid = GridNamelist(; mountain_height = 0.0E+1, mountain_half_width = 1.0E+4)
output = OutputNamelist(;
    output_variables = (:w,),
    output_file = output_file,
    output_steps = true,
)
sponge = SpongeNamelist(; use_sponge = false)
turbulence = TurbulenceNamelist(; turbulence_scheme = TKEScheme())

integrate(Namelists(; atmosphere, domain, grid, output, sponge, turbulence, setting))
