function compute_meanflow_effect!(state::State)
    if steady_state || single_column
        calc_integrals!(state, AbstractWKBMode())
    else
        calc_integrals!(state, MultiColumn())
    end

    set_boundaries!(state, BoundaryGWIntegrals())

    compute_gw_tendencies!(state)

    set_boundaries!(state, BoundaryGWTendencies())

    smooth_tendencies!(state)

    set_boundaries!(state, BoundaryGWTendencies())

    compute_gw_forcing!(state)

    set_boundaries!(state, BoundaryGWForces())

    compute_gw_heating!(state)

    return
end
