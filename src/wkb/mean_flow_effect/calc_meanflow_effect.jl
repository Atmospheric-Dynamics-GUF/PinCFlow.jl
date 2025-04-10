function calc_meanflow_effect!(state::State)
    if steady_state || single_column
        calc_integrals!(state, AbstractWKBMode())
    else
        calc_integrals!(state, Transient())
    end

    set_boundaries!(state, BoundaryGWIntegrals())

    calc_gw_tendencies!(state)

    set_boundaries!(state, BoundaryGWTendencies())

    smooth!(state)

    set_boundaries!(state, BoundaryGWTendencies())

    calc_gwmomforce!(state)

    set_boundaries!(state, BoundaryGWForces())

    calc_gwh!(state)

    return
end
