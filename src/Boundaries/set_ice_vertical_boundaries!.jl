"""
```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
```

Enforce vertical boundary conditions for tracers by dispatching to the appropriate method.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::TracerOn,
)
```

Enforce vertical boundary conditions for tracer predictands.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::TracerOn,
)
```

Enforce vertical boundary conditions for tracer reconstructions.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracer_setup::TracerOn,
)
```

Enforce vertical boundary conditions for vertical tracer fluxes.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    tracer_setup::TracerOn,
)
```

Enforce vertical boundary conditions for tracer WKB quantities by dispatching to the appropriate method.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
```

Enforce vertical boundary conditions for tracer WKB integrals.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
```

Enforce vertical boundary conditions for tracer WKB tendencies.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
"""
function set_ice_vertical_boundaries! end

function set_ice_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
    (; ice_setup) = state.namelists.ice
    set_ice_vertical_boundaries!(state, variables, ice_setup)
    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    ice_setup::NoIce,
)
    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    ice_setup::IceOn,
)
    (; namelists, domain) = state
    (; icepredictands) = state.ice

    for field in fieldnames(IcePredictands)
        set_vertical_boundaries_of_field!(
            getfield(icepredictands, field),
            namelists,
            domain,
            -,
        )
    end

    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    ice_setup::IceOn,
)
    (; namelists, domain) = state
    (; icereconstructions) = state.ice

    for field in fieldnames(IceReconstructions)
        set_vertical_boundaries_of_field!(
            getfield(icereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    ice_setup::IceOn,
)
    (; z_size) = state.namelists.domain
    (; nz, ko, k0, k1) = state.domain
    (; icefluxes) = state.ice
    @ivy if ko == 0
        for field in fieldnames(IceFluxes)
            getfield(icefluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
    end

    @ivy if ko + nz == z_size
        for field in fieldnames(IceFluxes)
            getfield(icefluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    ice_setup::IceOn,
)
    (; wkb_mode) = state.namelists.wkb
    set_ice_vertical_boundaries!(state, variables, wkb_mode)
    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
    println("Setting ice vertical boundaries for WKB integrals not finished")
    exit(1)
    (; namelists, domain) = state
    (; chiq0) = state.ice.iceforcings

    for field in (:uchi, :vchi, :wchi)
        set_vertical_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain,
            +;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
    println("Setting ice vertical boundaries for WKB integrals not finished")
    exit(1)
    (; namelists, domain) = state
    (; chiq0) = state.tracer.tracerforcings

    for field in (:dchidt,)
        set_vertical_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain,
            +,
        )
    end

    return
end
