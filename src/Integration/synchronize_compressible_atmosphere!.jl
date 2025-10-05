"""
```julia
synchronize_compressible_atmosphere!(state::State, predictands::Predictands)
```

Synchronize `state.atmosphere.pbar` with `predictands.p` if the atmosphere is compressible by dispatching to the appropriate method.

```julia
synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
)
```

Return in non-compressible modes.

```julia
synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
    model::Compressible,
)
```

Synchronize `state.atmosphere.pbar` with `predictands.p`.

In compressible mode, ``P`` is time-dependent. In the update of ``P``, only `state.variables.predictands.p` is changed, so that the old values of ``P`` are retained in the respective field of `state.atmosphere`. When these are no longer needed, this method is used to update the field accordingly.

# Arguments

  - `state`: Model state.

  - `predictands`: Predictands to use for the synchronization of the mass-weighted potential temperature.

  - `model`: Dynamic equations.
"""
function synchronize_compressible_atmosphere! end

function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
)
    (; model) = state.namelists.setting
    synchronize_compressible_atmosphere!(state, predictands, model)
    return
end

function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
)
    return
end

function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
    model::Compressible,
)
    (; pbar) = state.atmosphere
    (; p) = predictands

    pbar .= p

    return
end
