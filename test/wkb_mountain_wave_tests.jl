domain = DomainNamelist(;
    sizex = 40,
    sizey = 40,
    sizez = 40,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (0.0, 4.0E+5),
    ly_dim = (0.0, 4.0E+5),
    lz_dim = (0.0, 2.0E+4),
    nprocx = 1,
    nprocy = 1,
)

output = OutputNamelist(;
    output_variables = (),
    prepare_restart = true,
    restart = false,
    iin = -1,
    output_steps = true,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 3.6E+3,
    maxtime = 3.6E+3,
    fancy_namelists = true,
    input_file = "./pincflow_input.h5",
    output_file = "./pincflow_output.h5",
)

setting = SettingNamelist(;
    model = PseudoIncompressible(),
    testcase = WKBMountainWave(),
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
    mountainheight_dim = 1.5E+2,
    mountainwidth_dim = 5.0E+3,
    mountain_case = 13,
    height_factor = 2.0E+0,
    width_factor = 1.0E+1,
    spectral_modes = 1,
    stretch_exponent = 1.0E+0,
)

sponge = SpongeNamelist(;
    spongelayer = true,
    sponge_uv = false,
    spongeheight = 1.0E-1,
    spongealphaz_dim = 1.79E-2,
    spongealphaz_fac = 1.0E+0,
    unifiedsponge = true,
    lateralsponge = true,
    spongetype = ExponentialSponge(),
    spongeorder = 1,
    cosmosteps = 1,
    relax_to_mean = false,
    relaxation_period = 0.0E+0,
    relaxation_amplitude = 0.0E+0,
)

wkb = WKBNamelist(;
    xrmin_dim = 0.0E+0,
    xrmax_dim = 4.0E+5,
    yrmin_dim = 0.0,
    yrmax_dim = 4.0E+5,
    zrmin_dim = 0.0,
    zrmax_dim = 2.0E+4,
    nrxl = 1,
    nryl = 1,
    nrzl = 1,
    nrk_init = 1,
    nrl_init = 1,
    nrm_init = 1,
    nray_fac = 4,
    fac_dk_init = 1.0E-1,
    fac_dl_init = 1.0E-1,
    fac_dm_init = 1.0E-1,
    branchr = -1,
    merge_mode = ConstantWaveAction(),
    nsmth_wkb = 2,
    lsmth_wkb = true,
    sm_filter = Shapiro(),
    zmin_wkb_dim = 0.0E+0,
    lsaturation = true,
    alpha_sat = 1.0E+0,
    wkb_mode = MultiColumn(),
    blocking = false,
    nwm = 1,
    launch_algorithm = Clip(),
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
    wkb = wkb,
)

integrate(namelists)

# data = h5open("pincflow_output.h5")
# reference = h5open("wkb_mountain_wave_tests.h5")

# for key in keys(reference)
#     @test all(isapprox.(data[key], reference[key]))
# end

# close(data)
# close(reference)

# rm("pincflow_output.h5")
