"""
```julia
set_meridional_boundaries!(
    state::State,
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
)::Nothing
```

Enforce meridional boundary conditions for predictands or reconstructions by dispatching to the appropriate method.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:Boussinesq},
)::Nothing
```

Enforce meridional boundary conditions for predictands in Boussinesq mode.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:PseudoIncompressible},
)::Nothing
```

Enforce meridional boundary conditions for predictands in pseudo-incompressible mode.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:Compressible},
)::Nothing
```

Enforce meridional boundary conditions for predictands in compressible mode.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Val{:Boussinesq},
)::Nothing
```

Enforce meridional boundary conditions for reconstructions in Boussinesq mode.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Union{Val{:PseudoIncompressible}, Val{:Compressible}},
)::Nothing
```

Enforce meridional boundary conditions for reconstructions in non-Boussinesq modes.

```julia
set_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)::Nothing
```

Enforce meridional boundary conditions for WKB variables by dispatching to the appropriate method.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}},
)::Nothing
```

Enforce meridional boundary conditions for WKB integrals needed in `:SingleColumn` and `:SteadyState` configurations.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Val{:MultiColumn},
)::Nothing
```

Enforce meridional boundary conditions for WKB integrals needed in `:MultiColumn` configurations.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}},
)::Nothing
```

Enforce meridional boundary conditions for WKB tendencies needed in `:SingleColumn` and `:SteadyState` configurations.

```julia
set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Val{:MultiColumn},
)::Nothing
```

Enforce meridional boundary conditions for WKB tendencies needed in `:MultiColumn` configurations.

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
)::Nothing
    (; model) = state.namelists.atmosphere
    @dispatch_model set_meridional_boundaries!(state, variables, Val(model))
    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:Boussinesq},
)::Nothing
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in (:rhop, :u, :v, :w, :pip)
        set_meridional_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:PseudoIncompressible},
)::Nothing
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in (:rho, :rhop, :u, :v, :w, :pip)
        set_meridional_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Val{:Compressible},
)::Nothing
    (; namelists, domain) = state
    (; predictands) = state.variables

    for field in fieldnames(Predictands)
        set_meridional_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Val{:Boussinesq},
)::Nothing
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in (:rhoptilde, :utilde, :vtilde, :wtilde)
        set_meridional_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Union{Val{:PseudoIncompressible}, Val{:Compressible}},
)::Nothing
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in fieldnames(Reconstructions)
        set_meridional_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)::Nothing
    (; wkb_mode) = state.namelists.wkb
    @dispatch_wkb_mode set_meridional_boundaries!(
        state,
        variables,
        Val(wkb_mode),
    )
    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}},
)::Nothing
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e)
        set_meridional_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Val{:MultiColumn},
)::Nothing
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :utheta, :vtheta, :e)
        set_meridional_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}},
)::Nothing
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt)
        set_meridional_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    nothing
end

function set_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Val{:MultiColumn},
)::Nothing
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dthetadt)
        set_meridional_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    nothing
end
