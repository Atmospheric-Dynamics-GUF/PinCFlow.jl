"""
```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
```

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
```

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
```

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
```
"""
function set_turbulence_zonal_boundaries! end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; turbulencepredictands) = state.turbulence

    for field in fieldnames(TurbulencePredictands)
        set_zonal_boundaries_of_field!(
            getfield(turbulencepredictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; turbulencereconstructions) = state.turbulence

    for field in fieldnames(TurbulenceReconstructions)
        set_zonal_boundaries_of_field!(
            getfield(turbulencereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end
