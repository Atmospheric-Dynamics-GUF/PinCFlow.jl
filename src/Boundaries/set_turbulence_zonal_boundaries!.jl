"""
```julia
set_turbulence_zonal_boundaries!(state::State, variables::AbstractBoundaryVariables)
```

Enforce zonal boundary conditions for turbulences by dispatching to a turbulence-configuration-specific method.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::NoTurbulence,
)
```

Return for configurations without turbulence transport.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::AbstractTurbulence,
)
```

Enforce zonal boundary conditions for turbulences.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::NoTurbulence,
)
```

Return for configurations without turbulence transport.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::AbstractTurbulence,
)
```

Enforce zonal boundary conditions for reconstructions of turbulences.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
    turbulence_scheme::NoTurbulence,
)
```

Return for configurations without turbulence transport.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
    turbulence_scheme::AbstractTurbulence,
)
```

Enforce zonal boundary conditions for turbulence-gravity-wave-integral fields.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    turbulence_scheme::NoTurbulence,
)
```

Return for configurations without turbulence transport.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    turbulence_scheme::AbstractTurbulence,
)
```

Enforce zonal boundary conditions for turbulence-gravity-wave-tendency fields.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulence_scheme`: General turbulence-transport configuration.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_turbulence_zonal_boundaries! end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
    (; turbulence_scheme) = state.namelists.turbulence
    set_turbulence_zonal_boundaries!(state, variables, turbulence_scheme)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::NoTurbulence,
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::AbstractTurbulence,
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
    turbulence_scheme::NoTurbulence,
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::AbstractTurbulence,
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

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
    turbulence_scheme::NoTurbulence,
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
    turbulence_scheme::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; chiq0) = state.turbulence.turbulenceforcings

    for field in (:uchi, :vchi, :wchi)
        set_zonal_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    turbulence_scheme::NoTurbulence,
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    turbulence_scheme::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; chiq0) = state.turbulence.turbulenceforcings

    for field in (:dchidt,)
        set_zonal_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end
