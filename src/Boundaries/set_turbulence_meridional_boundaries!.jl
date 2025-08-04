"""
```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
```
"""
function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
    return
end

"""
```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
```
"""
function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; turbulencepredictands) = state.turbulence

    for field in fieldnames(TurbulencePredictands)
        set_meridional_boundaries_of_field!(
            getfield(turbulencepredictands, field),
            namelists,
            domain,
        )
    end

    return
end

"""
```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
```
"""
function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
    return
end

"""
```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
```
"""
function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; turbulencereconstructions) = state.turbulence

    for field in fieldnames(TurbulenceReconstructions)
        set_meridional_boundaries_of_field!(
            getfield(turbulencereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end
