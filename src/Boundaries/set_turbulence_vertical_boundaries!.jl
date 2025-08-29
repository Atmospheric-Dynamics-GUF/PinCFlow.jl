"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
```

Return for configurations without turbulence physics.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
```

Enforce vertical boundary conditions for prognostic turbulence variables.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
```

Return for configurations without turbulence physics.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
```

Enforce vertical boundary conditions for reconstructions of turbulence variables.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulencesetup::NoTurbulence,
)
```

Return for configurations without turbulence physics.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulencesetup::AbstractTurbulence,
)
```

Set the vertical fluxes of turbulence variables at the vertical boundaries to zero.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulencesetup`: General turbulence-physics configuration.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
"""
function set_turbulence_vertical_boundaries! end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; turbulencepredictands) = state.turbulence

    for field in fieldnames(TurbulencePredictands)
        set_vertical_boundaries_of_field!(
            getfield(turbulencepredictands, field),
            namelists,
            domain,
            -,
        )
    end

    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; turbulencereconstructions) = state.turbulence

    for field in fieldnames(TurbulenceReconstructions)
        set_vertical_boundaries_of_field!(
            getfield(turbulencereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulencesetup::NoTurbulence,
)
    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulencesetup::AbstractTurbulence,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; turbulencefluxes) = state.turbulence

    if ko == 0
        for field in fieldnames(TurbulenceFluxes)
            getfield(turbulencefluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
    end

    if ko + nzz == sizezz
        for field in fieldnames(TurbulenceFluxes)
            getfield(turbulencefluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end
