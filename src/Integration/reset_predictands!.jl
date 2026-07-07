"""
```julia
reset_predictands!(
    state::State,
    predictands::Predictands,
    tracerpredictands::TracerPredictands,
)
```

Reset fields in `state` to those in `predictands` and `tracerpredictands` by dispatching to specific methods.

```julia
reset_predictands!(state::State, tracerpredictands::TracerPredictands)
```

Reset fields in `state.tracer.tracerpredictands` to those in `tracerpredictands`.

```julia
reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Union{Val{:Boussinesq}, Val{:PseudoIncompressible}},
)
```

Reset the density, density fluctuations and wind components in `state.variables.predictands` to those in `predictands`.

```julia
reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:Compressible},
)
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
)
    (; model) = state.namelists.atmosphere

    @dispatch_model reset_predictands!(state, predictands, Val(model))
    reset_predictands!(state, tracerpredictands)

    return
end

function reset_predictands!(state::State, tracerpredictands::TracerPredictands)
    (; chi) = state.tracer.tracerpredictands

    @share chi = tracerpredictands.chi

    return
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:Boussinesq},
)
    (; rhop, u, v, w) = state.variables.predictands

    @share rhop = predictands.rhop
    @share u = predictands.u
    @share v = predictands.v
    @share w = predictands.w

    return
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:PseudoIncompressible},
)
    (; rho, rhop, u, v, w) = state.variables.predictands

    @share rho = predictands.rho
    @share rhop = predictands.rhop
    @share u = predictands.u
    @share v = predictands.v
    @share w = predictands.w

    return
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:Compressible},
)
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands

    @share rho = predictands.rho
    @share rhop = predictands.rhop
    @share u = predictands.u
    @share v = predictands.v
    @share w = predictands.w
    @share pip = predictands.pip
    @share p = predictands.p

    return
end
