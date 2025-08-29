"""
```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::NoIce,
)
```

Return for configurations without ice physics.

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::AbstractIce,
)
```

Enforce vertical boundary conditions for prognostic ice variables.

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::NoIce,
)
```

Return for configurations without ice physics.

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::AbstractIce,
)
```

Enforce vertical boundary conditions for reconstructions of ice variables.

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    icesetup::NoIce,
)
```

Return for configurations without ice physics.

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    icesetup::AbstractIce,
)
```

Set the vertical fluxes of ice variables at the vertical boundaries to zero.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `icesetup`: General ice-physics configuration.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
"""
function set_ice_vertical_boundaries! end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::NoIce,
)
    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::AbstractIce,
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
    icesetup::NoIce,
)
    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::AbstractIce,
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
    icesetup::NoIce,
)
    return
end

function set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    icesetup::AbstractIce,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; icefluxes) = state.ice

    if ko == 0
        for field in fieldnames(IceFluxes)
            getfield(icefluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
    end

    if ko + nzz == sizezz
        for field in fieldnames(IceFluxes)
            getfield(icefluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end
