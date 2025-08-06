"""
```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::NoIce,
)
```

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::AbstractIce,
)
```

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::NoIce,
)
```

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::AbstractIce,
)
```

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    icesetup::NoIce,
)
```

```julia
set_ice_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    icesetup::AbstractIce,
)
```
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
    (; zboundaries) = namelists.setting
    (; icepredictands) = state.ice

    for field in fieldnames(IcePredictands)
        set_vertical_boundaries_of_field!(
            getfield(icepredictands, field),
            namelists,
            domain,
            zboundaries,
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
    (; zboundaries) = namelists.setting
    (; icereconstructions) = state.ice

    for field in fieldnames(IceReconstructions)
        set_vertical_boundaries_of_field!(
            getfield(icereconstructions, field),
            namelists,
            domain,
            zboundaries,
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
