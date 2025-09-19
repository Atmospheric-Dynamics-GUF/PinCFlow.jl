"""
```julia
explicit_integration!(
    state::State,
    p0::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    side::LHS,
)
```

Integrate the left-hand sides of the prognostic equations with a Runge-Kutta time step.

At each Runge-Kutta stage, the prognostic variables are first reconstructed and their advective and diffusive fluxes are calculated. Subsequently, each variable is updated with its integrated left-hand side, followed immediately by an implicit Euler step (the size of which is the fractional time step at the current Runge-Kutta stage) that accounts for the Rayleigh-damping imposed by the unified sponge. After the Runge-Kutta loop, if the atmosphere is compressible, the Exner-pressure is updated with a full implicit Euler step to account for it being impacted by the Rayleigh damping of the mass-weighted potential temperature.

```julia
explicit_integration!(
    state::State,
    p0::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    side::RHS,
)
```

Perform an explicit Euler step on the right-hand sides of the prognostic equations and, if the atmosphere is compressible, update the Exner-pressure such that it is synchronized with the mass-weighted potential temperature.

# Arguments

  - `state`: Model state.

  - `p0`: The predictands that are used to compute the transporting velocities in the computation of the fluxes.

  - `dtstage`: Fractional time step.

  - `time`: Simulation time.

  - `side`: Side of the equations.
"""
function explicit_integration! end

function explicit_integration!(
    state::State,
    p0::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    side::LHS,
)
    (; nstages, stepfrac) = state.time
    (; tracersetup) = state.namelists.tracer
    (; icesetup) = state.namelists.ice
    (; turbulencesetup) = state.namelists.turbulence

    for rkstage in 1:nstages
        reconstruct!(state)
        set_boundaries!(state, BoundaryReconstructions())

        compute_fluxes!(state, p0)

        set_boundaries!(state, BoundaryFluxes())

        save_backups!(state, :rho)

        update!(state, dtstage, rkstage, Rho())
        apply_unified_sponge!(state, stepfrac[rkstage] * dtstage, time, Rho())

        update!(state, dtstage, rkstage, RhoP(), side)
        apply_unified_sponge!(state, stepfrac[rkstage] * dtstage, time, RhoP())

        update!(state, dtstage, rkstage, P())
        apply_unified_sponge!(state, stepfrac[rkstage] * dtstage, time, P())

        update!(state, dtstage, rkstage, tracersetup)
        apply_unified_sponge!(
            state,
            stepfrac[rkstage] * dtstage,
            time,
            tracersetup,
        )

        update!(state, dtstage, rkstage, icesetup)
        apply_unified_sponge!(
            state,
            stepfrac[rkstage] * dtstage,
            time,
            icesetup,
        )

        #update!(state, dtstage, rkstage, turbulencesetup)
        #apply_unified_sponge!(
            #state,
            #stepfrac[rkstage] * dtstage,
            #time,
            #turbulencesetup,
        #)

        set_boundaries!(state, BoundaryPredictands())

        update!(state, dtstage, rkstage, U(), side)
        update!(state, dtstage, rkstage, V(), side)
        update!(state, dtstage, rkstage, W(), side)
        apply_unified_sponge!(state, stepfrac[rkstage] * dtstage, time, U())
        apply_unified_sponge!(state, stepfrac[rkstage] * dtstage, time, V())
        apply_unified_sponge!(state, stepfrac[rkstage] * dtstage, time, W())

        set_boundaries!(state, BoundaryPredictands())
    end

    synchronize_compressible_atmosphere!(state, state.variables.predictands)

    apply_unified_sponge!(state, dtstage, time, PiP())

    return
end

function explicit_integration!(
    state::State,
    p0::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    side::RHS,
)
    synchronize_compressible_atmosphere!(state, state.variables.predictands)

    modify_compressible_wind!(state, *)

    set_boundaries!(state, BoundaryPredictands())

    save_backups!(state, :rhop, :u, :v, :w)

    update!(state, dtstage, RhoP(), RHS(), Explicit())

    update!(state, dtstage, U(), RHS(), Explicit())
    update!(state, dtstage, V(), RHS(), Explicit())
    update!(state, dtstage, W(), RHS(), Explicit())

    update!(state, dtstage, PiP())

    modify_compressible_wind!(state, /)

    set_boundaries!(state, BoundaryPredictands())

    return
end
