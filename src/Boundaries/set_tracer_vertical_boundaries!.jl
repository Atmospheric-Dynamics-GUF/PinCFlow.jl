"""
```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::AbstractTracer,
)
```

Enforce vertical boundary conditions for tracers.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::AbstractTracer,
)
```

Enforce vertical boundary conditions for reconstructions of tracers.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracersetup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracersetup::AbstractTracer,
)
```

Set the vertical tracer fluxes at the vertical boundaries to zero.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
  - `tracersetup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
"""
function set_tracer_vertical_boundaries! end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; tracerpredictands) = state.tracer

    for field in fieldnames(TracerPredictands)
        set_vertical_boundaries_of_field!(
            getfield(tracerpredictands, field),
            namelists,
            domain,
            zboundaries,
            -,
        )
    end

    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; tracerreconstructions) = state.tracer

    for field in fieldnames(TracerReconstructions)
        set_vertical_boundaries_of_field!(
            getfield(tracerreconstructions, field),
            namelists,
            domain,
            zboundaries,
        )
    end

    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracersetup::AbstractTracer,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; tracerfluxes) = state.tracer

    if ko == 0
        for field in fieldnames(TracerFluxes)
            getfield(tracerfluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
    end

    if ko + nzz == sizezz
        for field in fieldnames(TracerFluxes)
            getfield(tracerfluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end
