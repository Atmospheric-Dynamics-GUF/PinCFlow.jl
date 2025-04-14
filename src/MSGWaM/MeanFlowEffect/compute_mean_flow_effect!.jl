function compute_mean_flow_effect!(state::State)
    (; testcase) = state.namelists.setting
    compute_mean_flow_effect!(state, testcase)
    return
end

function compute_mean_flow_effect!(state::State, testcase::AbstractTestCase)
    return
end

function compute_mean_flow_effect!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb

    compute_integrals!(state, wkb_mode)

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
