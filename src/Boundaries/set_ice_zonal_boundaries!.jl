"""
```julia
set_ice_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::NoIce,
)
```

Return for configurations without ice physics.

```julia
set_ice_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::AbstractIce,
)
```

Enforce zonal boundary conditions for prognostic ice variables.

```julia
set_ice_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::NoIce,
)
```

Return for configurations without ice physics.

```julia
set_ice_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::AbstractIce,
)
```

Enforce zonal boundary conditions for reconstructions of ice variables.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `icesetup`: General ice-physics configuration.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_ice_zonal_boundaries! end

function set_ice_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::NoIce,
)
    return
end

function set_ice_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::AbstractIce,
)
    (; namelists, domain) = state
    (; icepredictands) = state.ice

    for field in fieldnames(IcePredictands)
        set_zonal_boundaries_of_field!(
            getfield(icepredictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_ice_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::NoIce,
)
    return
end

function set_ice_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::AbstractIce,
)
    (; namelists, domain) = state
    (; icereconstructions) = state.ice

    for field in fieldnames(IceReconstructions)
        set_zonal_boundaries_of_field!(
            getfield(icereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end
