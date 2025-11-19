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

h0 = 100.0
l0 = 1000.0

lx = 20000.0
ly = 20000.0
lz = 20000.0
dxr = lx / 2
dyr = ly / 2
dzr = lz / 2
alpharmax = 0.0179

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
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
<<<<<<< HEAD
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
<<<<<<< HEAD
>>>>>>> f0d2b4e
    npx = 8,
    npy = 8,
>>>>>>> 2aee3f7
=======
=======
    lx,
    ly,
    lz,
>>>>>>> cf395edbf2
    npx,
    npy,
    npz,
>>>>>>> afc93468
)
grid = GridNamelist(;
    resolved_topography = (x, y) -> h0 / (1 + (x^2 + y^2) / l0^2),
)
output =
    OutputNamelist(; output_variables = (:w,), output_file = "mountain_wave.h5")
sponge = SpongeNamelist(;
<<<<<<< HEAD
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
=======
    lhs_sponge = (x, y, z, t, dt) -> begin
        alpharx =
            abs(x) >= (lx - dxr) / 2 ?
            sin(pi * (abs(x) - (lx - dxr) / 2) / dxr)^2 : 0.0
        alphary =
            abs(y) >= (ly - dyr) / 2 ?
            sin(pi * (abs(y) - (ly - dyr) / 2) / dyr)^2 : 0.0
        alpharz =
            z >= lz - dzr ? sin(pi / 2 * (z - (lz - dzr)) / dzr)^2 : 0.0
        return alpharmax * (alpharx + alphary + alpharz) / 3
    end,
    relaxed_u = (x, y, z, t, dt) -> 10.0,
>>>>>>> cf395edbf2
)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))
