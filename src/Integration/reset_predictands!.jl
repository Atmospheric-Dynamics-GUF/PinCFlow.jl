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
    (; nxx, nyy, nzz) = state.domain
    (; rhop, u, v, w) = state.variables.predictands

    @share for k in 1:nzz, j in 1:nyy, i in 1:nxx
        rhop[i, j, k] = predictands.rhop[i, j, k]
        u[i, j, k] = predictands.u[i, j, k]
        v[i, j, k] = predictands.v[i, j, k]
        w[i, j, k] = predictands.w[i, j, k]
    end

    return
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:PseudoIncompressible},
)
    (; nxx, nyy, nzz) = state.domain
    (; rho, rhop, u, v, w) = state.variables.predictands

    @share for k in 1:nzz, j in 1:nyy, i in 1:nxx
        rho[i, j, k] = predictands.rho[i, j, k]
        rhop[i, j, k] = predictands.rhop[i, j, k]
        u[i, j, k] = predictands.u[i, j, k]
        v[i, j, k] = predictands.v[i, j, k]
        w[i, j, k] = predictands.w[i, j, k]
    end

    return
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Val{:Compressible},
)
    (; nxx, nyy, nzz) = state.domain
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands

    @share for k in 1:nzz, j in 1:nyy, i in 1:nxx
        rho[i, j, k] = predictands.rho[i, j, k]
        rhop[i, j, k] = predictands.rhop[i, j, k]
        u[i, j, k] = predictands.u[i, j, k]
        v[i, j, k] = predictands.v[i, j, k]
        w[i, j, k] = predictands.w[i, j, k]
        pip[i, j, k] = predictands.pip[i, j, k]
        p[i, j, k] = predictands.p[i, j, k]
    end

    return
end
