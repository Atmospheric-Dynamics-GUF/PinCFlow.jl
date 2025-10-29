# examples/submit/wkb_mountain_wave.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

h0 = 150.0
l0 = 5000.0
rl = 10
rh = 2

lx = 400000.0
ly = 400000.0
lz = 20000.0
dxr = lx / 20
dyr = ly / 20
dzr = lz / 10
alpharmax = 0.0179

atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0,
    initial_u = (x, y, z) -> 10.0,
)
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 40,
    lx,
    ly,
    lz,
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
    lhs_sponge = (x, y, z, t, dt) ->
        alpharmax / 3 * (
            exp((abs(x) - lx / 2) / dxr) +
            exp((abs(y) - ly / 2) / dyr) +
            exp((z - lz) / dzr)
        ),
    relaxed_u = (x, y, z, t, dt) -> 10.0,
)
wkb = WKBNamelist(; wkb_mode = MultiColumn())

integrate(Namelists(; atmosphere, domain, grid, output, sponge, wkb))
