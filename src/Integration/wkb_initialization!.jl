function wkb_initialization! end

function wkb_initialization!(state::State)
    (; tracer_setup, next_order_impact) = state.namelists.tracer

    initialize_rays!(state)

    if tracer_setup === :TracerOn && next_order_impact
        compute_gw_integrals!(state)
        set_boundaries!(state, BoundaryWKBIntegrals())
        smoothing!(state, Integrals())
    end
    return
end