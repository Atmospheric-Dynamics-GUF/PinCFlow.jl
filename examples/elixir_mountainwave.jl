using PinCFlow_dev

grid_list = setup_grid_list(nx = 150,  # Cells in x
                            ny = 1,    # Cells in y
                            nz = 50,  # Cells in z
                            nbx = 3,    # Halo/ghost cells in x
                            nby = 3,    # Halo/ghost cells in y
                            nbz = 3,    # Ghost cells in z
                            lx_dim = [0.0E+0, 6.0E+4],  # min, max of x
                            ly_dim = [0.0E+0, 4.0E+4],  # min, max of y
                            lz_dim = [0.0E+0, 2.0E+4],  # min, max of z
                            nprocx = 1,  # Processors in x
                            nprocy = 1)

output_list = setup_output_list(atmvarOut = "w",  # Atmospheric output variable
                                prepare_restart = false, # Save everything needed for a restart
                                restart = false, # Restart the model from previous simulation
                                iIn = -1,    # Restart at time step iIn
                                runName = "mountainwave", # Run name for netCDF file
                                outputType = "time", # 'timeStep' or 'time'
                                nOutput = 1,    # Output every nOutput time steps for outputType = 'timeStep'
                                maxIter = 1,    # Stop after maxIter time steps for outputType = 'timeStep'
                                outputTimeDiff = 3.6E+3, # Output every outputTimeDiff seconds for outputType = 'time'
                                maxTime = 3.6E+3, # Stop after maxTime seconds for outputType = 'time'
                                fancy_namelists = true)

debugging_list = setup_debugging_list(dtMin_dim = 1.0E-5)

test_case_list = setup_test_case_list(testCase = "mountainwave")

model_list = setup_model(model = "pseudo_incompressible")

solver_list = setup_solver_list(cfl = 5.0E-1, # CFL number
                                dtMax_dim = 6.0E+1, # Maximum time step
                                tStepChoice = "cfl", # 'fix' or 'cfl'
                                limiterType1 = "MCVariant")

poisson_solver_list = setup_poisson_solver_list(tolPoisson = 1.0E-8, # Abort criterion
                                                maxIterPoisson = 5000,   # Maximum iterations
                                                preconditioner = "yes", # 'no' or 'yes'
                                                dtau = 4.0E+0, # Time parameter for preconditioner
                                                maxIterADI = 2,      # Preconditioner iterations
                                                initialCleaning = true,   # Enforce initial non-divergence
                                                correctMomentum = true,   # Correct momentum to fulfill divergence constraint
                                                tolcrit = "abs")

atmosphere_list = setup_atmosphere_list(specifyReynolds = false, # Use inverse Reynolds number
                                        ReInv = 0.0E+0, # Inverse Reynolds number
                                        mu_viscous_dim = 0.0E+0, # Kinematic viscosity
                                        background = "isothermal", # 'isothermal'
                                        Temp0_dim = 3.0E+2, # Background temperature
                                        press0_dim = 1.0E+5, # Ground pressure
                                        backgroundFlow_dim = (1.0E+1, 0.0E+0, 0.0E+0), # Initial wind components
                                        f_Coriolis_dim = 0.0E+0, # Coriolis frequency
                                        corset = "constant")

topography_list = setup_topography_list(mountainHeight_dim = 4.0E+2, # Maximum height
                                        mountainWidth_dim = 1.0E+3, # Half width
                                        mountain_case = 3,      # Predefined topography
                                        range_factor = 1.0E+1, # Ratio between large and small scales
                                        spectral_modes = 1,      # Number of spectral modes
                                        envelope_reduction = 0.0E+0, # Relative reduction of the envelope (0 to 1)
                                        stretch_exponent = 1.0E+0)

boundary_list = setup_boundary_list(spongeLayer = true,  # General sponge layer switch
                                    sponge_uv = false, # Sponge layer for horizontal wind
                                    spongeHeight = 5.0E-1, # Relative height of lower sponge layer edge
                                    spongeAlphaZ_dim = 1.0E-2, # Maximum relaxation rate
                                    spongeAlphaZ_fac = 1.0E+0, # Sponge layer factor
                                    unifiedSponge = true,  # Unified sponge for both time schemes
                                    lateralSponge = true,  # Lateral sponge
                                    spongeType = "sinusoidal", # Sponge layer profile
                                    spongeOrder = 1,     # Order of polynomial sponge
                                    cosmoSteps = 1,     # Relative strength of COSMO sponge
                                    relax_to_mean = true,  # Relax the wind to its horizontal mean
                                    relaxation_period = 0.0E+0, # Period of an oscillation
                                    relaxation_amplitude = 0.0E+0, # Relative amplitude of oscillation
                                    xBoundary = "periodic", # Boundary conditions in x
                                    yBoundary = "periodic", # Boundary conditions in y
                                    zBoundary = "solid_wall")

semi = setup_semidiscretization(;
                                atmosphere_list, grid_list, output_list, debugging_list,
                                test_case_list, model_list, solver_list,
                                poisson_solver_list, topography_list, boundary_list);

dt = 30.0

pincflow(semi, dt);
