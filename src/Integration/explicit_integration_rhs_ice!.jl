function explicit_integration_rhs_ice!(
    state::State,
    dtstage::AbstractFloat,
)
    (; nstages) = state.time
    (; tracersetup) = state.namelists.tracer
    (; icesetup, dt_ice) = state.namelists.ice
    (; tref) = state.constants

    n_step_ice = ceil(Int, dtstage * tref / dt_ice)
    dtt_ice = dtstage / n_step_ice   

    for ii in 1:n_step_ice
        for rkstage in 1:nstages

        #reconstruct!(state)
        #set_boundaries!(state, BoundaryReconstructions())

        #compute_fluxes!(state, p0)

        #set_boundaries!(state, BoundaryFluxes())

        #save_backups!(state, :rho)
        ComputeSourceIce!(state) 
        update!(state, dtt_ice, rkstage, IceUpdatePhy())
        
        #apply_unified_sponge!(
        #    state,
        #    stepfrac[rkstage] * dtstage,
        #    time,
        #    tracersetup,
        #)

        set_boundaries!(state, BoundaryPredictands())

        end
    end

    return
end

