"""
```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
```

Return for configurations without turbulence physics.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
```

Enforce zonal boundary conditions for prognostic turbulence variables.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
```

Return for configurations without turbulence physics.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
```

Enforce zonal boundary conditions for reconstructions of turbulence variables.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulencesetup`: General turbulence-physics configuration.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
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
