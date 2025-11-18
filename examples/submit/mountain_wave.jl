# examples/submit/mountain_wave.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

<<<<<<< HEAD
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
=======
npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1
>>>>>>> afc93468

atmosphere = AtmosphereNamelist(;
    initial_wind = (1.0E+1, 0.0E+0, 0.0E+0),
    coriolis_frequency = 0.0E+0,
)
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
<<<<<<< HEAD
>>>>>>> f0d2b4e
    npx = 8,
    npy = 8,
>>>>>>> 2aee3f7
=======
    npx,
    npy,
    npz,
>>>>>>> afc93468
)
grid = GridNamelist(;
    resolved_topography = (x, y) -> 100 / (1 + (x^2 + y^2) / 1000^2),
)
output =
    OutputNamelist(; output_variables = (:w,), output_file = "mountain_wave.h5")
sponge = SpongeNamelist(;
<<<<<<< HEAD
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
=======
>>>>>>> afc93468
    alpharmax = 1.79E-2,
    lateral_sponge = true,
    sponge_type = SinusoidalSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
>>>>>>> 2aee3f7
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
