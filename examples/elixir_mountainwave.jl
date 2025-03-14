using PinCFlow

domain_list = DomainParameters(sizex=5, sizez=5)
output_list = OutputParameters()
debugging_list = DebugParameters()
test_case_list = TestCaseParameters()
model_list = ModelParameters()
discretization_list = DiscretizationParameters()
poisson_solver_list = PoissonSolverParameters()
atmosphere_list = AtmosphereParameters()
topography_list = TopographyParameters()
boundary_list = BoundaryParameters()

pars = Parameters(domain_list, output_list, debugging_list,
    test_case_list, model_list, discretization_list,
    poisson_solver_list, atmosphere_list, topography_list, boundary_list)

m = Model(pars)

dt = 30.0
pincflow(m, dt)
