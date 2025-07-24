"""
    set_boundaries!(state, variables::BoundaryPredictands)

Set all boundary conditions (zonal, meridional, vertical) for predictand fields.
"""
function set_boundaries!(state::State, variables::BoundaryPredictands)
    (; zboundaries) = state.namelists.setting
    (; tracersetup) = state.namelists.tracer 

    divide_tracer_by_rho!(state, tracersetup)

    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)

    multiply_tracer_by_rho!(state, tracersetup)
    return
end

function divide_tracer_by_rho!(state::State, tracersetup::NoTracer)
    return
end

function divide_tracer_by_rho!(state::State, tracersetup::AbstractTracer)

    (; chi) = state.tracer.tracerpredictands
    (; rho) = state.variables.predictands
    (; rhostrattfc) = state.atmosphere

    chi .= chi ./ (rho .+ rhostrattfc)

    return 
end

function multiply_tracer_by_rho!(state::State, tracersetup::NoTracer)
    return 
end

function multiply_tracer_by_rho!(state::State, tracersetup::AbstractTracer)

    (; chi) = state.tracer.tracerpredictands
    (; rho) = state.variables.predictands
    (; rhostrattfc) = state.atmosphere

    chi .= chi .* (rho .+ rhostrattfc)

    return 
end

"""
    set_boundaries!(state, variables::BoundaryReconstructions)

Set all boundary conditions for reconstruction fields.
"""
function set_boundaries!(state::State, variables::BoundaryReconstructions)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

"""
    set_boundaries!(state, variables::BoundaryFluxes)

Set vertical boundary conditions for flux fields (no horizontal boundaries needed).
"""
function set_boundaries!(state::State, variables::BoundaryFluxes)
    (; zboundaries) = state.namelists.setting
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

"""
    set_boundaries!(state, variables::BoundaryGWIntegrals)

Set all boundary conditions for gravity wave integral fields.
"""
function set_boundaries!(state::State, variables::BoundaryGWIntegrals)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

"""
    set_boundaries!(state, variables::BoundaryGWTendencies)

Set all boundary conditions for gravity wave tendency fields.
"""
function set_boundaries!(state::State, variables::BoundaryGWTendencies)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end
