"""
```julia
set_meridional_boundaries!(
    state::State,
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
)
```

Enforce meridional boundary conditions for predictands or reconstructions by dispatching to the appropriate method.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Boussinesq,
)
```

Enforce meridional boundary conditions for predictands in Boussinesq mode.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::PseudoIncompressible,
)
```

Enforce meridional boundary conditions for predictands in pseudo-incompressible mode.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
```

Enforce meridional boundary conditions for predictands in compressible mode.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Boussinesq,
)
```

Enforce meridional boundary conditions for reconstructions in Boussinesq mode.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Union{PseudoIncompressible, Compressible},
)
```

Enforce meridional boundary conditions for reconstructions in non-Boussinesq modes.

```julia
set_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
```

Enforce meridional boundary conditions for WKB variables by dispatching to the appropriate method.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn},
)
```

Enforce meridional boundary conditions for WKB integrals needed in `SingleColumn` and `SteadyState` configurations.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::MultiColumn,
)
```

Enforce meridional boundary conditions for WKB integrals needed in `MultiColumn` configurations.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{SteadyState, SingleColumn},
)
```

Enforce meridional boundary conditions for WKB tendencies needed in `SingleColumn` and `SteadyState` configurations.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::MultiColumn,
)
```

Enforce meridional boundary conditions for WKB tendencies needed in `MultiColumn` configurations.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function set_meridional_boundaries! end

function set_meridional_boundaries!(
    state::State,
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
)
    (; model) = state.namelists.atmosphere
    set_meridional_boundaries!(state, variables, model)
    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Boussinesq,
)
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in (:rhop, :u, :v, :w, :pip)
        set_meridional_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::PseudoIncompressible,
)
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in (:rho, :rhop, :u, :v, :w, :pip)
        set_meridional_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in fieldnames(Predictands)
        set_meridional_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Boussinesq,
)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in (:rhoptilde, :utilde, :vtilde, :wtilde)
        set_meridional_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Union{PseudoIncompressible, Compressible},
)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in fieldnames(Reconstructions)
        set_meridional_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
    (; wkb_mode) = state.namelists.wkb
    set_meridional_boundaries!(state, variables, wkb_mode)
    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn},
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e, :sterm, :bterm)
        set_meridional_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :utheta, :vtheta, :e, :sterm, :bterm)
        set_meridional_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{SteadyState, SingleColumn},
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dtkedt)
        set_meridional_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dthetadt, :dtkedt)
        set_meridional_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryDiffusionCoefficients,
)
    (; namelists, domain) = state
    (; turbulencediffusioncoefficients) = state.turbulence

    for field in (:kh, :km, :kek)
        set_meridional_boundaries_of_field!(
            getfield(turbulencediffusioncoefficients, field),
            namelists,
            domain,
        )
    end

    return
end
