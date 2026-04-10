"""
```julia
set_zonal_boundaries!(
    state::State,
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
)
```

Enforce zonal boundary conditions for predictands or reconstructions by dispatching to the appropriate method.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:Boussinesq},
)
```

Enforce zonal boundary conditions for predictands in Boussinesq mode.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:PseudoIncompressible},
)
```

Enforce zonal boundary conditions for predictands in pseudo-incompressible mode.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:Compressible},
)
```

Enforce zonal boundary conditions for predictands in compressible mode.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Val{:Boussinesq},
)
```

Enforce zonal boundary conditions for reconstructionss in Boussinesq mode.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Union{Val{:PseudoIncompressible}, Val{:Compressible}},
)
```

Enforce zonal boundary conditions for reconstructions in non-Boussinesq modes.

```julia
set_zonal_boundaries!(state::State, variables::AbstractBoundaryWKBVariables)
```

Enforce zonal boundary conditions for WKB variables by dispatching to the appropriate method.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}},
)
```

Enforce zonal boundary conditions for WKB integrals needed in `:SingleColumn` and `:SteadyState` configurations.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Val{:MultiColumn},
)
```

Enforce zonal boundary conditions for WKB integrals needed in `:MultiColumn` configurations.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}},
)
```

Enforce zonal boundary conditions for WKB tendencies needed in `:SingleColumn` and `:SteadyState` configurations.

```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Val{:MultiColumn},
)
```

Enforce zonal boundary conditions for WKB tendencies needed in `:MultiColumn` configurations.

```julia
set_zonal_boundaries!(state::State, variables::BoundaryDiffusionCoefficients)
```

Enforce zonal boundary conditions for eddy diffusion coefficients.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_zonal_boundaries! end

function set_zonal_boundaries!(
    state::State,
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
)
    (; model) = state.namelists.atmosphere
    @dispatch_model set_zonal_boundaries!(state, variables, Val(model))
    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:Boussinesq},
)
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in (:rhop, :u, :v, :w, :pip)
        set_zonal_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:PseudoIncompressible},
)
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in (:rho, :rhop, :u, :v, :w, :pip)
        set_zonal_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:Compressible},
)
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in fieldnames(Predictands)
        set_zonal_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Val{:Boussinesq},
)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in (:rhoptilde, :utilde, :vtilde, :wtilde)
        set_zonal_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Union{Val{:PseudoIncompressible}, Val{:Compressible}},
)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in fieldnames(Reconstructions)
        set_zonal_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
    (; wkb_mode) = state.namelists.wkb
    @dispatch_wkb_mode set_zonal_boundaries!(state, variables, Val(wkb_mode))
    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}},
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e, :uhat, :vhat, :what, :bhat, :pihat)
        set_zonal_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Val{:MultiColumn},
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :utheta, :vtheta, :e)
        set_zonal_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    for field in (:uhat, :vhat, :what, :bhat, :pihat)
        set_zonal_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}},
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt)
        set_zonal_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Val{:MultiColumn},
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dthetadt)
        set_zonal_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end

function set_zonal_boundaries!(
    state::State,
    variables::BoundaryDiffusionCoefficients,
)
    (; namelists, domain) = state
    (; turbulencediffusioncoefficients) = state.turbulence
    for field in (:kh, :km, :kek)
        set_zonal_boundaries_of_field!(
            getfield(turbulencediffusioncoefficients, field),
            namelists,
            domain,
        )
    end

    return
end
