function set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)
    (; namelists, domain) = state
    (; rho, rhop, u, v, w, pip) = state.variables.predictands

    set_meridional_boundaries_of_field!(rho, namelists, domain)
    set_meridional_boundaries_of_field!(rhop, namelists, domain)
    set_meridional_boundaries_of_field!(u, namelists, domain)
    set_meridional_boundaries_of_field!(v, namelists, domain)
    set_meridional_boundaries_of_field!(w, namelists, domain)
    set_meridional_boundaries_of_field!(pip, namelists, domain)

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
    (; namelists, domain) = state
    (; rhotilde, rhoptilde, utilde, vtilde, wtilde) =
        state.variables.reconstructions

    set_meridional_boundaries_of_field!(rhotilde, namelists, domain)
    set_meridional_boundaries_of_field!(rhoptilde, namelists, domain)
    set_meridional_boundaries_of_field!(utilde, namelists, domain)
    set_meridional_boundaries_of_field!(vtilde, namelists, domain)
    set_meridional_boundaries_of_field!(wtilde, namelists, domain)

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
)
    (; wkb_mode) = state.namelists.wkb
    set_meridional_boundaries!(state, variables, wkb_mode)
    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn},
)
    (; namelists, domain) = state
    (; uw, vw, e) = state.wkb.integrals

    set_meridional_boundaries_of_reduced_field!(uw, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(vw, namelists, domain)

    # This one might be unnecessary.
    set_meridional_boundaries_of_reduced_field!(e, namelists, domain)

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; uu, uv, uw, vv, vw, etx, ety, e, utheta, vtheta) = state.wkb.integrals

    set_meridional_boundaries_of_reduced_field!(uu, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(uv, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(uw, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(vv, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(vw, namelists, domain)

    # These might be unnecessary.
    set_meridional_boundaries_of_reduced_field!(e, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(etx, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(ety, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(utheta, namelists, domain)
    set_meridional_boundaries_of_reduced_field!(vtheta, namelists, domain)

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
)
    (; wkb_mode) = state.namelists.wkb
    set_meridional_boundaries!(state, variables, wkb_mode)
    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::Union{SteadyState, SingleColumn},
)
    (; namelists, domain) = state
    (; dudt, dvdt) = state.wkb.tendencies

    set_meridional_boundaries_of_field!(dudt, namelists, domain)
    set_meridional_boundaries_of_field!(dvdt, namelists, domain)

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; dudt, dvdt, dthetadt) = state.wkb.tendencies

    set_meridional_boundaries_of_field!(dudt, namelists, domain)
    set_meridional_boundaries_of_field!(dvdt, namelists, domain)
    set_meridional_boundaries_of_field!(dthetadt, namelists, domain)

    return
end
