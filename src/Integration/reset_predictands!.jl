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

# added ice predictands (ice active vars) reset
function reset_predictands!(state::State, source_icepredictands::IcePredictands) 
    # Get the current active ice state descriptor
    active_ice_descriptor = get_IceActiveVars(state.ice.icepredictands)
    
    # Extract the tuple of active fields
    active_fields = ice_active_vars_tuple(active_ice_descriptor)
    
    # Dynamically reset only the active fields
    for field_name in active_fields
        target_array = getfield(state.ice.icepredictands, field_name)
        source_array = getfield(source_icepredictands, field_name)
        target_array .= source_array
    end
    
    return
end
