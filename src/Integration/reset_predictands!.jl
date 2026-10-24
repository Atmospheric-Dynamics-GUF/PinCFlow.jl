"""
```julia
reset_predictands!(
    state::State,
    predictands::Predictands,
    tracerpredictands::TracerPredictands,
)::Nothing
```

Reset fields in `state` to those in `predictands` and `tracerpredictands` by dispatching to specific methods.

```julia
reset_predictands!(state::State, tracerpredictands::TracerPredictands)::Nothing
```

Reset fields in `state.tracer.tracerpredictands` to those in `tracerpredictands`.

```julia
reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Union{Val{:Boussinesq}, Val{:PseudoIncompressible}},
)::Nothing
```

Reset the density, density fluctuations and wind components in `state.variables.predictands` to those in `predictands`.

```julia
reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:Compressible},
)::Nothing
```

Reset the density, density fluctuations, wind components, Exner pressure and mass-weighted potential temperature (i.e. all fields) in `state.variables.predictands` to those in `predictands`.

# Arguments

  - `state`: Model state.

  - `predictands`: Fields to reset to.

  - `tracerpredictands`: Tracer fields to reset to.

  - `model`: Dynamic equations.
"""
function reset_predictands! end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    tracerpredictands::TracerPredictands,
)::Nothing
    (; model) = state.namelists.atmosphere

    @dispatch_model reset_predictands!(state, predictands, Val(model))
    reset_predictands!(state, tracerpredictands)

    nothing
end

function reset_predictands!(
    state::State,
    tracerpredictands::TracerPredictands,
)::Nothing
    (; chi) = state.tracer.tracerpredictands

    chi .= tracerpredictands.chi

    nothing
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:Boussinesq},
)::Nothing
    (; rhop, u, v, w) = state.variables.predictands

    rhop .= predictands.rhop
    u .= predictands.u
    v .= predictands.v
    w .= predictands.w

    nothing
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:PseudoIncompressible},
)::Nothing
    (; rho, rhop, u, v, w) = state.variables.predictands

    rho .= predictands.rho
    rhop .= predictands.rhop
    u .= predictands.u
    v .= predictands.v
    w .= predictands.w

    nothing
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:Compressible},
)::Nothing
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands

    rho .= predictands.rho
    rhop .= predictands.rhop
    u .= predictands.u
    v .= predictands.v
    w .= predictands.w
    pip .= predictands.pip
    p .= predictands.p

    nothing
end
