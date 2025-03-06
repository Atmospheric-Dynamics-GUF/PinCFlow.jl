# TODO - Reorder them as in the original namelist file

function setup_atmosphere_list(; specifyReynolds, # Use inverse Reynolds number
                               ReInv, # Inverse Reynolds number
                               mu_viscous_dim, # Kinematic viscosity
                               background, # 'isothermal'
                               Temp0_dim, # Background temperature
                               press0_dim, # Ground pressure
                               backgroundFlow_dim, # Initial wind components
                               f_Coriolis_dim, # Coriolis frequency
                               corset)
    return (; specifyReynolds, ReInv, mu_viscous_dim, background, Temp0_dim, press0_dim,
            backgroundFlow_dim, f_Coriolis_dim, corset)
end

function setup_grid_list(; nx, ny, nz, nbx, nby, nbz, lx_dim, ly_dim, lz_dim, nprocx,
                         nprocy)
    return (; nx, ny, nz, nbx, nby, nbz, lx_dim, ly_dim, lz_dim, nprocx, nprocy)
end

function setup_output_list(; atmvarOut, prepare_restart, restart, iIn, runName, outputType,
                           nOutput, maxIter, outputTimeDiff, maxTime, fancy_namelists)
    return (; atmvarOut, prepare_restart, restart, iIn, runName, outputType, nOutput,
            maxIter, outputTimeDiff, maxTime, fancy_namelists)
end

function setup_debugging_list(; dtMin_dim)
    return (; dtMin_dim)
end

function setup_test_case_list(; testCase)
    return (; testCase)
end

function setup_model(; model)
    return (; model)
end

function setup_solver_list(; cfl, dtMax_dim, tStepChoice, limiterType1)
    return (; cfl, dtMax_dim, tStepChoice, limiterType1)
end

function setup_poisson_solver_list(; tolPoisson, maxIterPoisson, preconditioner, dtau,
                                   maxIterADI, initialCleaning, correctMomentum, tolcrit)
    return (; tolPoisson, maxIterPoisson, preconditioner, dtau, maxIterADI, initialCleaning,
            correctMomentum, tolcrit)
end

function setup_topography_list(; mountainHeight_dim, mountainWidth_dim, mountain_case,
                               range_factor,
                               spectral_modes, envelope_reduction, stretch_exponent)
    return (; mountainHeight_dim, mountainWidth_dim, mountain_case, range_factor,
            spectral_modes,
            envelope_reduction, stretch_exponent)
end

function setup_boundary_list(; spongeLayer, sponge_uv, spongeHeight, spongeAlphaZ_dim,
                             spongeAlphaZ_fac, unifiedSponge, lateralSponge, spongeType,
                             spongeOrder, cosmoSteps,
                             relax_to_mean, relaxation_period, relaxation_amplitude,
                             xBoundary, yBoundary, zBoundary)
    return (; spongeLayer, sponge_uv, spongeHeight, spongeAlphaZ_dim, spongeAlphaZ_fac,
            unifiedSponge, lateralSponge, spongeType, spongeOrder, cosmoSteps,
            relax_to_mean,
            relaxation_period, relaxation_amplitude, xBoundary, yBoundary, zBoundary)
end

function setup_semidiscretization(; atmosphere_list, grid_list, output_list, debugging_list,
                                  test_case_list, model_list, solver_list,
                                  poisson_solver_list, topography_list, boundary_list)
    (; specifyReynolds, ReInv, mu_viscous_dim, background, Temp0_dim, press0_dim,
    backgroundFlow_dim, f_Coriolis_dim, corset) = atmosphere_list
    (; nx, ny, nz, nbx, nby, nbz, lx_dim, ly_dim, lz_dim, nprocx, nprocy) = grid_list
    (; atmvarOut, prepare_restart, restart, iIn, runName, outputType, nOutput, maxIter, outputTimeDiff, maxTime, fancy_namelists) = output_list
    (; dtMin_dim) = debugging_list
    (; testCase) = test_case_list
    (; model) = model_list
    (; cfl, dtMax_dim, tStepChoice, limiterType1) = solver_list
    (; tolPoisson, maxIterPoisson, preconditioner, dtau, maxIterADI, initialCleaning,
    correctMomentum, tolcrit) = poisson_solver_list
    (; mountainHeight_dim, mountainWidth_dim, mountain_case, range_factor, spectral_modes,
    envelope_reduction, stretch_exponent) = topography_list
    (; spongeLayer, sponge_uv, spongeHeight, spongeAlphaZ_dim, spongeAlphaZ_fac,
    unifiedSponge, lateralSponge, spongeType, spongeOrder, cosmoSteps, relax_to_mean,
    relaxation_period, relaxation_amplitude, xBoundary, yBoundary, zBoundary) = boundary_list

    semi = initialize_values(nx, ny, nz, nbx, nby, nbz, lx_dim[1], lx_dim[2], ly_dim[1],
                             ly_dim[2], lz_dim[1], lz_dim[2]; model, maxIterPoisson,
                             tolcrit,
                             tolPoisson, maxTime, tStepChoice, dtMax_dim, cfl, Temp0_dim,
                             specifyReynolds, mu_viscous_dim, background,
                             press0_dim, backgroundFlow_dim, f_Coriolis_dim, corset, ReInv,
                             stretch_exponent, spongeLayer, sponge_uv)

    return semi
end
