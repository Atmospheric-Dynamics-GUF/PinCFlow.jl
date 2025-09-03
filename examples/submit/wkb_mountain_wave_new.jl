# examples/submit/wkb_mountain_wave.jl

using PinCFlow

if length(ARGS) == 0
    output_file = "./pincflow_output.h5"
elseif length(ARGS) == 1
    output_file = ARGS[1] * "/pincflow_output.h5"
else
    error("Too many arguments to the script!")
end

atmosphere = AtmosphereNamelist(; backgroundflow_dim = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
    sizex = 5,
    sizey = 5,
    sizez = 5,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (-2.0E+5, 2.0E+5),
    ly_dim = (-2.0E+5, 2.0E+5),
    lz_dim = (0.0E+0, 2.0E+4),
    npx = 1,
    npy = 1,
)

grid = GridNamelist(;
    mountainheight_dim = 1.5E+2,
    mountainwidth_dim = 5.0E+3,
    mountain_case = 13,
    height_factor = 2.0E+0,
    width_factor = 1.0E+1,
)

output = OutputNamelist(; output_variables = (:w,), output_file = output_file)

setting =
    SettingNamelist(; model = Compressible(), testcase = WKBMountainWave())

sponge = SpongeNamelist(;
    spongelayer = true,
    spongeheight = 1.0E-1,
    spongealphaz_dim = 1.79E-2,
    unifiedsponge = true,
    lateralsponge = true,
    spongetype = ExponentialSponge(),
    relax_to_mean = false,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)

tracer =
    TracerNamelist(; tracersetup = LinearTracer(), leading_order_impact = true)

wkb = WKBNamelist(;
    xrmin_dim = -2.0E+5,
    xrmax_dim = 2.0E+5,
    yrmin_dim = -2.0E+5,
    yrmax_dim = 2.0E+5,
    zrmin_dim = 0.0,
    zrmax_dim = 2.0E+4,
)

integrate(
    Namelists(; atmosphere, domain, grid, output, setting, sponge, wkb, tracer),
)
