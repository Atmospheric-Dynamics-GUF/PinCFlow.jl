# examples/submit/wkb_mountain_wave.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

h0 = 150
l0 = 5000
rl = 10
rh = 2

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0E+0,
    initial_u = (x, y, z) -> 1.0E+1,
)
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 40,
    lx = 4.0E+5,
    ly = 4.0E+5,
    lz = 2.0E+4,
    npx,
    npy,
    npz,
)
grid = GridNamelist(;
    resolved_topography = (x, y) ->
        x^2 + y^2 <= (rl * l0)^2 ?
        h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) * rh / (rh + 1) : 0.0,
    unresolved_topography = (alpha, x, y) ->
        x^2 + y^2 <= (rl * l0)^2 ?
        (
            pi / l0,
            0.0,
            h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) / (rh + 1),
        ) : (0.0, 0.0, 0.0),
)
output = OutputNamelist(;
    output_variables = (:w,),
    output_file = "wkb_mountain_wave.h5",
)
sponge = SpongeNamelist(;
    sponge_extent = 1.0E-1,
    alpharmax = 1.79E-2,
    lateral_sponge = true,
    sponge_type = ExponentialSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)
wkb = WKBNamelist(; wkb_mode = MultiColumn())

integrate(Namelists(; atmosphere, domain, grid, output, sponge, wkb))
