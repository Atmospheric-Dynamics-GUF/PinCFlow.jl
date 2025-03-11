using PinCFlow_dev

domain_list = DomainParameters(sizex=10)
output_list = OutputParameters()
debugging_list = DebugParameters()
test_case_list = TestCaseParameters()
model_list = ModelParameters()
discretization_list = DiscretizationParameters()
poisson_solver_list = PoissonSolverParameters(tolpoisson=1.0E-8, # Abort criterion
    maxiterpoisson=5000,   # Maximum iterations
    preconditioner="yes", # 'no' or 'yes'
    dtau=4.0E+0, # Time parameter for preconditioner
    maxiteradi=2,      # Preconditioner iterations
    initalcleaning=true,
    correctmomentum=true,   # Correct momentum to fulfill divergence constraint
    tolcrit="abs")

atmosphere_list = AtmosphereParameters(specifyreynolds=false, # Use inverse Reynolds number
    reinv=0.0E+0, # Inverse Reynolds number
    mu_viscous_dim=0.0E+0, # Kinematic viscosity
    background="isothermal", # 'isothermal'
    temp0_dim=3.0E+2, # Background temperature
    press0_dim=1.0E+5, # Ground pressure
    backgroundflow_dim=[1.0E+1, 0.0E+0, 0.0E+0], # Initial wind components
    f_coriolis_dim=0.0E+0, # Coriolis frequency
    corset="constant")

topography_list = TopographyParameters(
    mountainheight_dim=4.0E+2, # Maximum height
    mountainwidth_dim=1.0e+3, # half width
    mountain_case=3,      # predefined topography
    range_factor=1.0e+1, # ratio between large and small scales
    spectral_modes=1,      # number of spectral modes
    envelope_reduction=0.0e+0, # relative reduction of the envelope (0 to 1)
    stretch_exponent=1.0e+0)

boundary_list = BoundaryParameters(
    spongelayer=true,  # general sponge layer switch
    sponge_uv=false, # sponge layer for horizontal wind
    spongeheight=5.0e-1, # relative height of lower sponge layer edge
    spongealphaz_dim=1.0e-2, # maximum relaxation rate
    spongealphaz_fac=1.0e+0, # sponge layer factor
    unifiedsponge=true,  # unified sponge for both time schemes
    lateralsponge=true,  # lateral sponge
    spongetype="sinusodial", # sponge layer profile
    spongeorder=1,     # order of polynomial sponge
    cosmosteps=1,     # relative strength of cosmo sponge
    relax_to_mean=true,  # relax the wind to its horizontal mean
    relaxation_period=0.0e+0, # period of an oscillation
    relaxation_amplitude=0.0e+0, # relative amplitude of oscillation
    xboundary="periodic", # boundary conditions in x
    yboundary="periodic", # boundary conditions in y
    zboundary="solid_wall")

pars = Parameters(domain_list, output_list, debugging_list,
    test_case_list, model_list, discretization_list,
    poisson_solver_list, atmosphere_list, topography_list, boundary_list)

m = Model(pars)

dt = 30.0
pincflow(m, dt)
