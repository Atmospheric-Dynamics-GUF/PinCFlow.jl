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
    model::Union{Boussinesq, PseudoIncompressible},
)
```

Reset the density, density fluctuations and wind components in `state.variables.predictands` to those in `predictands`.

```julia
reset_predictands!(state::State, predictands::Predictands, model::Compressible)
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
    icepredictands::IcePredictands,
)
    (; model) = state.namelists.atmosphere

    reset_predictands!(state, predictands, model)
    reset_predictands!(state, tracerpredictands)
    reset_predictands!(state, icepredictands) # added ice predictands reset

    return
end

function reset_predictands!(state::State, tracerpredictands::TracerPredictands)
    (; chi) = state.tracer.tracerpredictands

    chi .= tracerpredictands.chi

    return
end

function reset_predictands!(state::State, icepredictands::IcePredictands) # added ice predictands reset
    (; n, q, qv) = state.ice.icepredictands

    n .= icepredictands.n
    q .= icepredictands.q
    qv .= icepredictands.qv

    return
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Union{Boussinesq, PseudoIncompressible},
)
    (; rho, rhop, u, v, w) = state.variables.predictands

    rho .= predictands.rho
    rhop .= predictands.rhop
    u .= predictands.u
    v .= predictands.v
    w .= predictands.w

    return
end

function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::Compressible,
)
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands

    rho .= predictands.rho
    rhop .= predictands.rhop
    u .= predictands.u
    v .= predictands.v
    w .= predictands.w
    pip .= predictands.pip
    p .= predictands.p

    return
end
