include("../../../src/PinCFlow.jl")

using .PinCFlow

output_file =
    "./examples/submit/local/pincflow_output_mountain.h5"

domain = DomainNamelist(;
    sizex = 40,
    sizey = 40,
    sizez = 40,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (0.0, 2.0E+4),
    ly_dim = (0.0, 2.0E+4),
    lz_dim = (0.0, 2.0E+4),
    npx = 1,
    npy = 1,
)

output = OutputNamelist(;
    output_variables = (:w, :rhop, ),
    prepare_restart = false,
    restart = false,
    iin = -1,
    output_steps = false,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 1.0,
    maxtime = 1.0,
    fancy_namelists = true,
    input_file = "./pincflow_input.h5",
    output_file = output_file,
)

setting = SettingNamelist(;
    model = Boussinesq(),
    testcase = MountainWave(),
    zboundaries = SolidWallBoundaries(),
)

discretization = DiscretizationNamelist(;
    cfl = 5.0E-1,
    dtmin_dim = 1.0E-5,
    dtmax_dim = 6.0E+1,
    adaptive_time_step = true,
    limitertype = MCVariant(),
)

poisson = PoissonNamelist(;
    tolpoisson = 1.0E-8,
    maxiterpoisson = 1000,
    preconditioner = true,
    dtau = 4.0E+0,
    maxiteradi = 2,
    initialcleaning = true,
    relative_tolerance = false,
)

atmosphere = AtmosphereNamelist(;
    specifyreynolds = false,
    reinv = 0.0E+0,
    mu_viscous_dim = 0.0E+0,
    background = Isothermal(),
    temp0_dim = 3.0E+2,
    press0_dim = 1.0E+5,
    backgroundflow_dim = (1.0E+1, 0.0E+0, 0.0E+0),
    coriolis_frequency = 0.0E+0,
    coriolis_mode = FPlane(),
)

grid = GridNamelist(;
    mountainheight_dim = 1.0E+2,
    mountainwidth_dim = 1.0E+3,
    mountain_case = 4,
    height_factor = 1.0E+0,
    width_factor = 1.0E+0,
    spectral_modes = 1,
    stretch_exponent = 1.0E+0,
)

sponge = SpongeNamelist(;
    spongelayer = true,
    sponge_uv = false,
    spongeheight = 5.0E-1,
    spongealphaz_dim = 1.79E-2,
    spongealphaz_fac = 1.0E+0,
    unifiedsponge = true,
    lateralsponge = true,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
    spongetype = SinusoidalSponge(),
    spongeorder = 1,
    cosmosteps = 1,
    relax_to_mean = false,
    perturbation_period = 0.0E+0,
    perturbation_amplitude = 0.0E+0,
)

tracer = TracerNamelist(;
    tracersetup = LinearTracer(),
)

namelists = Namelists(;
    domain = domain,
    output = output,
    setting = setting,
    discretization = discretization,
    poisson = poisson,
    atmosphere = atmosphere,
    grid = grid,
    sponge = sponge,
    tracer = tracer,
)

integrate(namelists)
