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

    compute_gw_integrals!(state, wkb_mode)

    set_boundaries!(state, BoundaryGWIntegrals())

    compute_gw_tendencies!(state)

    set_boundaries!(state, BoundaryGWTendencies())

    smooth_gw_tendencies!(state)

    set_boundaries!(state, BoundaryGWTendencies())

    return
end
