function set_zonal_boundaries!(state::State, variables::BoundaryPredictands)
    (; namelists, domain) = state
    (; predictands) = state.variables
    (; model) = namelists.setting

    for field in (:rho, :rhop, :u, :v, :w, :pip)
        set_zonal_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    set_compressible_zonal_boundaries!(state, model)

    return
end

function set_zonal_boundaries!(state::State, variables::BoundaryReconstructions)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in fieldnames(Reconstructions)
        set_zonal_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(state::State, variables::BoundaryGWIntegrals)
    (; wkb_mode) = state.namelists.wkb
    set_zonal_boundaries!(state, variables, wkb_mode)
    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e)
        set_zonal_boundaries_of_reduced_field!(
            getfield(integrals, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :etx, :ety, :utheta, :vtheta, :e)
        set_zonal_boundaries_of_reduced_field!(
            getfield(integrals, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(state::State, variables::BoundaryGWTendencies)
    (; wkb_mode) = state.namelists.wkb
    set_zonal_boundaries!(state, variables, wkb_mode)
    return
end

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
