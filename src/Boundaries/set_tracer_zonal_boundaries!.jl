"""
```julia
set_tracer_zonal_boundaries!(state::State, variables::AbstractBoundaryVariables)
```

Enforce zonal boundary conditions for tracers by dispatching to the appropriate method.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::Val{:NoTracer},
)
```

Return for configurations without tracer transport.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::Val{:TracerOn},
)
```

Enforce zonal boundary conditions for tracer predictands.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::Val{:TracerOn},
)
```

Enforce zonal boundary conditions for tracer reconstructions.

```julia 
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracer_setup::Val{:TracerOn},
)
```

Return for tracer fluxes, as these only need vertical boundary conditions enforced.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    tracer_setup::Val{:TracerOn},
)
```

Enforce zonal boundary conditions for tracer WKB integrals.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    tracer_setup::Val{:TracerOn},
)
```

Enforce zonal boundary conditions for tracer WKB tendencies.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_tracer_zonal_boundaries! end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
    (; tracer_setup) = state.namelists.tracer
    @dispatch_tracer_setup set_tracer_zonal_boundaries!(
        state,
        variables,
        Val(tracer_setup),
    )
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::Val{:NoTracer},
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::Val{:TracerOn},
)
    (; namelists, domain) = state
    (; tracerpredictands) = state.tracer

    for field in fieldnames(TracerPredictands)
        set_zonal_boundaries_of_field!(
            getfield(tracerpredictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::Val{:TracerOn},
)
    (; namelists, domain) = state
    (; tracerreconstructions) = state.tracer

    for field in fieldnames(TracerReconstructions)
        set_zonal_boundaries_of_field!(
            getfield(tracerreconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracer_setup::Val{:TracerOn},
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    tracer_setup::Val{:TracerOn},
)
    (; namelists, domain) = state
    (; tracerwkbintegrals) = state.tracer
    (; leading_order_impact) = namelists.tracer

    if leading_order_impact
        for field in (:uchi0, :vchi0, :wchi0)
            set_zonal_boundaries_of_field!(
                getfield(tracerwkbintegrals, field),
                namelists,
                domain;
                layers = (1, 1, 1),
            )
        end
    end

    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    tracer_setup::Val{:TracerOn},
)
    (; namelists, domain) = state
    (; dchidt0) = state.tracer.tracerwkbtendencies
    (; leading_order_impact) = namelists.tracer

    if leading_order_impact
        set_zonal_boundaries_of_field!(
            dchidt0,
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end
