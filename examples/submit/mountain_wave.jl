# examples/submit/mountain_wave.jl

using PinCFlow

<<<<<<< HEAD
if length(ARGS) == 0
	output_file = "./pincflow_output.h5"
=======
@ivy if length(ARGS) == 0
    output_file = "./pincflow_output.h5"
>>>>>>> 2aee3f7
elseif length(ARGS) == 1
	output_file = ARGS[1] * "/pincflow_output.h5"
else
	error("Too many arguments to the script!")
end

atmosphere = AtmosphereNamelist(; initial_wind = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
<<<<<<< HEAD
<<<<<<< HEAD
	sizex = 40,
	sizey = 40,
	sizez = 40,
	lx_dim = (-1.0E+4, 1.0E+4),
	ly_dim = (-1.0E+4, 1.0E+4),
	lz_dim = (0.0E+0, 2.0E+4),
	npx = 8,
	npy = 8,
=======
    sizex = 40,
    sizey = 40,
    sizez = 40,
    lx_dim = 2.0E+4,
    ly_dim = 2.0E+4,
    lz_dim = 2.0E+4,
=======
    x_size = 40,
    y_size = 40,
    z_size = 40,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
>>>>>>> f0d2b4e
    npx = 8,
    npy = 8,
>>>>>>> 2aee3f7
)
grid = GridNamelist(; mountain_case = 4)
output = OutputNamelist(; output_variables = (:w,), output_file = output_file)
sponge = SpongeNamelist(;
<<<<<<< HEAD
<<<<<<< HEAD
	spongelayer = true,
	spongealphaz_dim = 1.79E-2,
	unifiedsponge = true,
	lateralsponge = true,
	spongetype = SinusoidalSponge(),
	relax_to_mean = false,
	relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
=======
    spongelayer = true,
=======
    use_sponge = true,
>>>>>>> f0d2b4e
    alpharmax = 1.79E-2,
    betarmax = 0.0E+0,
    lateral_sponge = true,
    sponge_type = SinusoidalSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
>>>>>>> 2aee3f7
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
