"""
    set_zonal_boundaries!(state, variables::BoundaryPredictands)

Set zonal boundaries for all predictand fields (rho, rhop, u, v, w, pip) and handle
compressible model boundaries.
"""
function set_zonal_boundaries!(state::State, variables::BoundaryPredictands)
    (; namelists, domain) = state
    (; predictands) = state.variables
    (; model) = namelists.setting
    (; tracersetup) = namelists.tracer
    (; icesetup) = namelists.ice
    (; turbulencesetup) = namelists.turbulence

    for field in (:rho, :rhop, :u, :v, :w, :pip)
        set_zonal_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    set_compressible_zonal_boundaries!(state, model)

    set_tracer_zonal_boundaries!(state, variables, tracersetup)
    set_ice_zonal_boundaries!(state, variables, icesetup)
    set_turbulence_zonal_boundaries!(state, variables, turbulencesetup)

    return
end

"""
    set_zonal_boundaries!(state, variables::BoundaryReconstructions)

Set zonal boundaries for all reconstruction fields.
"""
function set_zonal_boundaries!(state::State, variables::BoundaryReconstructions)
    (; namelists, domain) = state
    (; reconstructions) = state.variables
    (; tracersetup) = namelists.tracer
    (; icesetup) = namelists.ice
    (; turbulencesetup) = namelists.turbulence

    for field in fieldnames(Reconstructions)
        set_zonal_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    set_tracer_zonal_boundaries!(state, variables, tracersetup)
    set_ice_zonal_boundaries!(state, variables, icesetup)
    set_turbulence_zonal_boundaries!(state, variables, turbulencesetup)

    return
end

"""
    set_zonal_boundaries!(state, variables::BoundaryGWIntegrals)

Set zonal boundaries for gravity wave integral fields. Dispatches based on WKB mode.
"""
function set_zonal_boundaries!(state::State, variables::BoundaryGWIntegrals)
    (; wkb_mode) = state.namelists.wkb
    set_zonal_boundaries!(state, variables, wkb_mode)
    return
end

"""
    set_zonal_boundaries!(state, variables::BoundaryGWIntegrals, wkb_mode::AbstractWKBMode)

Set zonal boundaries for basic GW integral fields (uw, vw, e) with minimal boundary layers.
"""
function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e)
        set_zonal_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

"""
    set_zonal_boundaries!(state, variables::BoundaryGWIntegrals, wkb_mode::MultiColumn)

Set zonal boundaries for extended GW integral fields in multi-column mode, including
cross-correlations (uu, uv, vv) and energy transport terms.
"""
function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :etx, :ety, :utheta, :vtheta, :e)
        set_zonal_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

"""
    set_zonal_boundaries!(state, variables::BoundaryGWTendencies)

Set zonal boundaries for GW tendency fields. Dispatches based on WKB mode.
"""
function set_zonal_boundaries!(state::State, variables::BoundaryGWTendencies)
    (; wkb_mode) = state.namelists.wkb
    set_zonal_boundaries!(state, variables, wkb_mode)
    return
end

"""
    set_zonal_boundaries!(state, variables::BoundaryGWTendencies, wkb_mode::AbstractWKBMode)

Set zonal boundaries for basic GW tendency fields (dudt, dvdt).
"""
function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt)
        set_zonal_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end

"""
    set_zonal_boundaries!(state, variables::BoundaryGWTendencies, wkb_mode::MultiColumn)

Set zonal boundaries for GW tendency fields in multi-column mode, including
temperature tendency (dthetadt).
"""
function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dthetadt)
        set_zonal_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end
