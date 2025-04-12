function set_zonal_boundaries!(state::State, variables::BoundaryPredictands)
    (; namelists, domain) = state
    (; rho, rhop, u, v, w, pip) = state.variables.predictands

    set_zonal_boundaries_of_field!(rho, namelists, domain)
    set_zonal_boundaries_of_field!(rhop, namelists, domain)
    set_zonal_boundaries_of_field!(u, namelists, domain)
    set_zonal_boundaries_of_field!(v, namelists, domain)
    set_zonal_boundaries_of_field!(w, namelists, domain)
    set_zonal_boundaries_of_field!(pip, namelists, domain)

    return
end

function set_zonal_boundaries!(state::State, variables::BoundaryReconstructions)
    (; namelists, domain) = state
    (; rhotilde, rhoptilde, utilde, vtilde, wtilde) =
        state.variables.reconstructions

    set_zonal_boundaries_of_field!(rhotilde, namelists, domain)
    set_zonal_boundaries_of_field!(rhoptilde, namelists, domain)
    set_zonal_boundaries_of_field!(utilde, namelists, domain)
    set_zonal_boundaries_of_field!(vtilde, namelists, domain)
    set_zonal_boundaries_of_field!(wtilde, namelists, domain)

    return
end

function set_zonal_boundaries!(state::State, variables::BoundaryGWIntegrals)
    (; namelists, domain) = state
    (; uu, uv, uw, vv, vw, etx, ety, e, utheta, vtheta) = state.wkb.integrals

    if steady_state || single_column # not sure if this is the right condition
        # all remaining fluxes are 0
        set_zonal_boundaries_of_reduced_field!(uw, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(vw, namelists, domain)
        # not necessary to set boundaries of e
        set_zonal_boundaries_of_reduced_field!(e, namelists, domain)
    else
        # transient mode
        set_zonal_boundaries_of_reduced_field!(uu, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(uv, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(uw, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(vv, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(vw, namelists, domain)
        # these are not necessary
        set_zonal_boundaries_of_reduced_field!(e, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(etx, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(ety, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(utheta, namelists, domain)
        set_zonal_boundaries_of_reduced_field!(vtheta, namelists, domain)
    end
    return
end

function set_zonal_boundaries!(state::State, variables::BoundaryGWTendencies)
    (; namelists, domain) = state
    (; dudt, dvdt, dthetadt) = state.wkb.integrals

    set_zonal_boundaries_of_reduced_field!(dudt, namelists, domain)
    set_zonal_boundaries_of_reduced_field!(dvdt, namelists, domain)
    if !steady_state || !single_column
        set_zonal_boundaries_of_reduced_field!(dthetadt, namelists, domain)
    end
    return
end

function set_zonal_boundaries!(state::State, variables::BoundaryGWForces)
    (; namelists, domain) = state
    (; u, v, w) = state.wkb.gwmomforce

    set_zonal_boundaries_of_reduced_field!(u, namelists, domain)
    set_zonal_boundaries_of_reduced_field!(v, namelists, domain)
    set_zonal_boundaries_of_reduced_field!(w, namelists, domain)
    return
end
