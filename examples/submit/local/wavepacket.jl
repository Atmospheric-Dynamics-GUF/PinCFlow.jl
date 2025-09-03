include("../../../src/PinCFlow.jl")

using .PinCFlow

output_file = "./examples/submit/local/pincflow_output.h5"

domain = DomainNamelist(;
    sizex = 32,
    sizey = 1,
    sizez = 100,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (0.0, 9000.0E+3),
    ly_dim = (0.0, 300.0E+3),
    lz_dim = (0.0, 100.0E+3),
    npx = 1,
    npy = 1,
)

output = OutputNamelist(;
    output_variables = (:w, :u),
    prepare_restart = false,
    restart = false,
    iin = -1,
    output_steps = false,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 9000.0,
    maxtime = 90000.0,
    fancy_namelists = true,
    input_file = "./pincflow_input.h5",
    output_file = output_file,
)

setting = SettingNamelist(;
    model = Compressible(),
    testcase = WKBWavePacket(),
    zboundaries = SolidWallBoundaries(),
)

discretization = DiscretizationNamelist(;
    cfl = 5.0E-1,
    dtmin_dim = 1.0E+0,
    dtmax_dim = 1.0E+2,
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
    backgroundflow_dim = (0.0E+0, 0.0E+0, 0.0E+0),
    coriolis_frequency = 1.E-4,
    coriolis_mode = FPlane(),
)

grid = GridNamelist(;
    mountainheight_dim = 0.0,
    mountainwidth_dim = 1.0E+3,
    mountain_case = 4,
    height_factor = 1.0E+0,
    width_factor = 1.0E+0,
    spectral_modes = 1,
    stretch_exponent = 1.0E+0,
)

sponge = SpongeNamelist(;
    spongelayer = false,
    sponge_uv = false,
    spongeheight = 5.0E-1,
    spongealphaz_dim = 1.79E-2,
    spongealphaz_fac = 1.0E+0,
    unifiedsponge = false,
    lateralsponge = false,
    relaxation_wind = (0.0E+0, 0.0E+0, 0.0E+0),
    spongetype = SinusoidalSponge(),
    spongeorder = 1,
    cosmosteps = 1,
    relax_to_mean = false,
    perturbation_period = 0.0E+0,
    perturbation_amplitude = 0.0E+0,
)

wkb = WKBNamelist(;
    xrmin_dim = 0.0E+0,
    xrmax_dim = 9000.0E+3,
    yrmin_dim = 0.0,
    yrmax_dim = 300.0E+3,
    zrmin_dim = 0.E+0,
    zrmax_dim = 100.E+3,
    nrxl = 2,
    nryl = 1,
    nrzl = 2,
    nrk_init = 1,
    nrl_init = 1,
    nrm_init = 1,
    nray_fac = 1,
    fac_dk_init = 1.0E-1,
    fac_dl_init = 1.0E-1,
    fac_dm_init = 1.0E-4,
    branchr = 1,
    merge_mode = ConstantWaveAction(),
    nsmth_wkb = 2,
    lsmth_wkb = true,
    sm_filter = Shapiro(),
    zmin_wkb_dim = 0.0E+0,
    lsaturation = false,
    alpha_sat = 1.0E+0,
    wkb_mode = MultiColumn(),
    blocking = false,
    nwm = 1,
)

tracer = TracerNamelist(;
    tracersetup = LinearTracer(),
    leading_order_impact = true,
    )

wavepacket = WavePacketNamelist(;
    wavepacketdim = 3,
    lambdax_dim = 0.0,
    lambday_dim = 300.E+3,
    lambdaz_dim = 1.E+3,
    x0_dim = 4500.E+3,
    y0_dim = 0.0,
    z0_dim = 10.E+3,
    sigmax_dim = 1500.E+3,
    sigmay_dim = 0.0,
    sigmaz_dim = 5.E+3,
    a0 = 0.5,
    branch = 1,
    u0_jet_dim = 3.E+1,
    sigmaz_jet_dim = 4.E+3,
    z0_jet_dim = 12.E+3,
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
    wavepacket = wavepacket,
    wkb = wkb,
    tracer = tracer,
)

# state = State(namelists)

integrate(namelists)
