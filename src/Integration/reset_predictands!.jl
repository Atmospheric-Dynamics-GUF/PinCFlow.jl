"""
```julia
reset_predictands!(state::State, predictands::Predictands)
```

Reset fields in `state.variables.predictands` to those in `predictands`.

This function dispatches to the appropriate model-specific methods.

# Arguments

  - `state`: Model state.
  - `predictands`: Fields to reset to.
"""
function reset_predictands!(state::State, predictands::Predictands)
    (; model) = state.namelists.setting
    reset_predictands!(state, predictands, model)
    return
end

"""
```julia
reset_predictands!(state::State, predictands::Predictands, model::AbstractModel)
```

Reset the density, density fluctuations and wind components in `state.variables.predictands` to those in `predictands`.

# Arguments

  - `state`: Model state.
  - `predictands`: Fields to reset to.
  - `model`: Dynamic equations.
"""
function reset_predictands!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
)
    (; rho, rhop, u, v, w) = state.variables.predictands

    rho .= predictands.rho
    rhop .= predictands.rhop
    u .= predictands.u
    v .= predictands.v
    w .= predictands.w

    return
end

"""
```julia
reset_predictands!(state::State, predictands::Predictands, model::Compressible)
```

Reset the density, density fluctuations, wind components, Exner pressure and mass-weighted potential temperature (i.e. all fields) in `state.variables.predictands` to those in `predictands`.

# Arguments

  - `state`: Model state.
  - `predictands`: Fields to reset to.
  - `model`: Dynamic equations.
"""
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
