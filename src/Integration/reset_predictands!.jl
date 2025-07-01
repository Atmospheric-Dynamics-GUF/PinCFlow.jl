"""
```julia
reset_predictands!(state::State, predictands::Predictands)
```

Reset current predictand fields to specified values based on model type.

This function dispatches to the appropriate model-specific implementation
for restoring predictand variables from a given state.

# Arguments

  - `state::State`: Simulation state to modify
  - `predictands::Predictands`: Source values to restore from
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

Reset predictand fields for general model types.

Restores the basic set of predictand variables (density, velocity components)
from the provided source state. Used for non-compressible model types.

# Arguments

  - `state::State`: Simulation state to modify
  - `predictands::Predictands`: Source values containing rho, rhop, u, v, w
  - `model::AbstractModel`: General model type (Boussinesq, PseudoIncompressible, etc.)
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

Reset predictand fields for compressible model.

Restores all predictand variables including pressure fields that are specific
to compressible flow models.

# Arguments

  - `state::State`: Simulation state to modify
  - `predictands::Predictands`: Source values containing all compressible variables
  - `model::Compressible`: Compressible model type

# Details

Restores: rho, rhop, u, v, w, pip (pressure perturbation), p (total pressure)
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
