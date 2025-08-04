"""
```julia
set_ice_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::NoIce,
)
```
"""
function set_ice_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::NoIce,
)
    return
end

"""
```julia
set_ice_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::AbstractIce,
)
```
"""
function set_ice_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    icesetup::AbstractIce,
)
    (; namelists, domain) = state
    (; icepredictands) = state.ice

    for field in fieldnames(IcePredictands)
        set_meridional_boundaries_of_field!(
            getfield(icepredictands, field),
            namelists,
            domain,
        )
    end

    return
end

"""
```julia
set_ice_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::NoIce,
)
```
"""
function set_ice_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::NoIce,
)
    return
end

"""
```julia
set_ice_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::AbstractIce,
)
```
"""
function set_ice_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    icesetup::AbstractIce,
)
    (; namelists, domain) = state
    (; icereconstructions) = state.ice

    for field in fieldnames(IceReconstructions)
        set_meridional_boundaries_of_field!(
            getfield(icereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end
