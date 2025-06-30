"""
    set_meridional_boundaries!(state::State, variables::BoundaryPredictands)

Enforce meridional boundary conditions for all predictand fields.
"""
function set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)
    (; namelists, domain) = state
    (; predictands) = state.variables
    (; model) = namelists.setting

    for field in (:rho, :rhop, :u, :v, :w, :pip)
        set_meridional_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    set_compressible_meridional_boundaries!(state, model)

    return
end

"""
    set_meridional_boundaries!(state::State, variables::BoundaryReconstructions)

Enforce meridional boundary conditions for all reconstruction fields.
"""
function set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in fieldnames(Reconstructions)
        set_meridional_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

"""
    set_meridional_boundaries!(state::State, variables::BoundaryGWIntegrals)

Enforce meridional boundary conditions for gravity-wave-integral fields, dispatching based on WKB mode.
"""
function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
)
    (; wkb_mode) = state.namelists.wkb
    set_meridional_boundaries!(state, variables, wkb_mode)
    return
end

"""
    set_meridional_boundaries!(state::State, variables::BoundaryGWIntegrals, wkb_mode::AbstractWKBMode)

Enforce meridional boundary conditions for gravity-wave-integral fields needed in `SingleColumn` and `SteadyState` configurations.
"""
function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e)
        set_meridional_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

"""
    set_meridional_boundaries!(state::State, variables::BoundaryGWIntegrals, wkb_mode::MultiColumn)

Enforce meridional boundary conditions for gravity-wave-integral fields needed in `MultiColumn` configurations.
"""
function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :etx, :ety, :utheta, :vtheta, :e)
        set_meridional_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

"""
    set_meridional_boundaries!(state::State, variables::BoundaryGWTendencies)

Enforce meridional boundary conditions for gravity-wave-tendency fields, dispatching based on WKB mode.
"""
function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
)
    (; wkb_mode) = state.namelists.wkb
    set_meridional_boundaries!(state, variables, wkb_mode)
    return
end

"""
    set_meridional_boundaries!(state::State, variables::BoundaryGWTendencies, wkb_mode::AbstractWKBMode)

Enforce meridional boundary conditions for gravity-wave-tendency fields needed in `SingleColumn` and `SteadyState` configurations.
"""
function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt)
        set_meridional_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end

"""
    set_meridional_boundaries!(state::State, variables::BoundaryGWTendencies, wkb_mode::MultiColumn)

Enforce meridional boundary conditions for gravity-wave-tendency fields needed in `MultiColumn` configurations.
"""
function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dthetadt)
        set_meridional_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end
