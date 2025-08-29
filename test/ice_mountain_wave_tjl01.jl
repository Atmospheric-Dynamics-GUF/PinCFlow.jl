#julia --project=. test/ice_mountain_wave.jl
#julia --project=. examples/visualization/fast_plot_compare.jl 
#using Revise

include("../src/PinCFlow.jl")

using .PinCFlow
#using HDF5

domain = DomainNamelist(;
    sizex = 40,
    sizey = 1,
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
    output_variables = (:w,:n,:qv,:q),
    prepare_restart = true,
    restart = false,
    iin = -1,
    output_steps = false,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 3.6E+3,
    maxtime = 3.6E+3,
    input_file = "./test/pincflow_input.h5",
    output_file = "./test/pincflow_output.h5",
)

setting = SettingNamelist(;
    model = PseudoIncompressible(),
    testcase = MountainWave(),
    zboundaries = SolidWallBoundaries(),
)

discretization = DiscretizationNamelist(;
    cfl = 5.0E-1,
    dtmin_dim = 1.0E-5,
    dtmax_dim = 1.0E+0,
    adaptive_time_step = true,
    limitertype = MCVariant(),
)

poisson = PoissonNamelist(;
    tolpoisson = 1.0E-8,
    maxiterpoisson = 1000,
    preconditioner = true,
    dtau = 4.0E+0,
    maxiteradi = 10,
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
    mountain_case = 3,
    height_factor = 1.0E+0,
    width_factor = 1.0E+0,
    spectral_modes = 1,
    stretch_exponent = 1.0E+0,
)

#sponge = SpongeNamelist(;
#    spongelayer = true,
#    sponge_uv = false,
#    spongeheight = 5.0E-1,
#    spongealphaz_dim = 1.79E-2,
#    spongealphaz_fac = 1.0E+0,
#    unifiedsponge = true,
#    lateralsponge = true,
#    spongetype = SinusoidalSponge(),
#    spongeorder = 1,
#    cosmosteps = 1,
#    relax_to_mean = false,
#    perturbation_period = 0.0E+0,
#    perturbation_amplitude = 0.0E+0,
#    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
#)
sponge = SpongeNamelist(;
    spongelayer = false,
    sponge_uv = false,
    spongeheight = 5.0E-1,
    spongealphaz_dim = 1.79E-2,
    spongealphaz_fac = 1.0E+0,
    unifiedsponge = false,
    lateralsponge = false,
    spongetype = SinusoidalSponge(),
    spongeorder = 1,
    cosmosteps = 1,
    relax_to_mean = true,
    perturbation_period = 0.0E+0,
    perturbation_amplitude = 0.0E+0,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)
ice = IceNamelist(
    icesetup = IceOn(),
    dt_ice = 1., 
    #ice_source = 1,
    #ice_predictands = IcePredictands(),
    #ice_auxiliaries = IceAuxiliaries(IcePredictands()),
    #ice_physics = IcePhysics(),
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
    ice=ice
)

integrate(namelists)

#rm("./test/pincflow_output.h5")174 
