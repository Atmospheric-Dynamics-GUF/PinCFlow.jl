#using Revise

include("../src/PinCFlow.jl")

using .PinCFlow
using HDF5

domain = DomainNamelist(;
    sizex = 8,
    sizey = 1,
    sizez = 8,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (0.0, 4.0E+4),
    ly_dim = (0.0, 1.0E+4),
    lz_dim = (0.0, 1.5E+4),
    npx = 1,    
    npy = 1,
    npz = 1,
)

output = OutputNamelist(;
    output_variables = (:rhop, :pip, :w, :u, :thetap, :n2, :rhobar, :thetabar, :n, :qv, :q, :iaux1, :iaux2, :iaux3, :wwp, :epp, :thp),
    prepare_restart = true,
    restart = false,
    iin = -1,
    output_steps = false,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 1.0, # 3.6E+1, #E+3
    maxtime = 12.0, #3.6E+1, #E+3
    input_file = "./test/pincflow_input.h5",
    output_file = "./test/pincflow_output.h5",
)

setting = SettingNamelist(;
    model = PseudoIncompressible(),
    #testcase = WKBMountainWave(),
    testcase = WKBMultipleWavePackets(),
    #testcase = MultipleWavePackets(),
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
    initialcleaning = false,
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
#    backgroundflow_dim = (1.0E+1, 0.0E+0, 0.0E+0), # mountain wave
    coriolis_frequency = 1.0E-4,
    coriolis_mode = FPlane(),
)

grid = GridNamelist(;
    #mountainheight_dim = 1.0E+3,
    mountainheight_dim = 0.0E+3,
    mountainwidth_dim = 1.0E+3,
    mountain_case = 3,
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
    spongetype = SinusoidalSponge(),
    spongeorder = 1,
    cosmosteps = 1,
    relax_to_mean = false,
    perturbation_period = 0.0E+0,
    perturbation_amplitude = 0.0E+0,
    relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)
ice = IceNamelist(;
    icesetup = IceOn(),
    dt_ice = 1.,
    nscx = 16,
    nscy = 1,
    nscz = 40,
    # nscx = 1,
    # nscy = 1,
    # nscz = 1,
    cloudcover = CloudCoverOn(),
#    large_scale_ice = true,
#    parameterized_nucleation = false,
 ) 
wkb = WKBNamelist(;
                    xrmin_dim = 0.0E+4,
                    xrmax_dim = 4.0E+4,
                    yrmin_dim = 0.0E+4,
                    yrmax_dim = 1.0E+4,
                    zrmin_dim = 0.0E+0,
                    zrmax_dim = 1.5E+4,
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
                    branchr = 1,
                    merge_mode = ConstantWaveAction(),
                    nsmth_wkb = 2,
                    #lsmth_wkb = true,
                    lsmth_wkb = false,
                    sm_filter = Shapiro(),
                    zmin_wkb_dim = 0.0,
                    lsaturation = false,
                    #lsaturation = true,
                    alpha_sat = 1.0E+0,
                    wkb_mode = MultiColumn(),
                    #blocking = true,
                    blocking = false,
                    long_threshold = 2.5E-1,
                    drag_coefficient = 1.0E+0,
                    nwm = 2,
                )

#multiwavepackets = MultiWavePacketNamelist(RandomWavePackets(), 2, domain, setting.testcase)
multiwavepackets = MultiWavePacketNamelist(;
    nwm = 2,
    wavepacketdim = [1, 1],
    lambdax_dim = [1.0E+4, 1.0E+4],
    lambday_dim = [0.0E+0   , 0.0E+0   ],
    lambdaz_dim = [-2.0E+3, 2.0E+3],
    x0_dim = [5.0E+3, 5.0E+3],
    y0_dim = [5.0E+3, 5.0E+3],
    z0_dim = [7.0E+3, 9.0E+3],
    sigmax_dim = [0.0E+3, 0.0E+3],
    sigmay_dim = [0.0E+3, 0.0E+3],
    sigmaz_dim = [6.0E+3, 6.0E+3],
    a0 = [0.12E+0, 0.12E+0],
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
    ice=ice,
    wkb=wkb,
#    wavepacket=wavepacket,
    multiwavepackets=multiwavepackets,
)

integrate(namelists)

#include("../examples/visualization/fast_plot_tjl04.jl")
