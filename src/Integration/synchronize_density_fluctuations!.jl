"""
```julia
synchronize_density_fluctuations!(state::State)
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho`.

This function dispatches to the appropriate model-specific methods.

# Arguments

  - `state`: Model state.
"""
function synchronize_density_fluctuations!(state::State)
    (; model) = state.namelists.setting
    synchronize_density_fluctuations!(state, model)
    return
end

"""
```julia
synchronize_density_fluctuations!(state::State, model::Boussinesq)
```

Return for Boussinesq mode.

In Boussinesq mode, density fluctuations don't require synchronization,
since the density is assumed constant except in the buoyancy equation.

# Arguments

  - `state`: Model state.
  - `model`: Dynamic equations.
"""
function synchronize_density_fluctuations!(state::State, model::Boussinesq)
    return
end

"""
```julia
synchronize_density_fluctuations!(state::State, model::PseudoIncompressible)
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho`.

In pseudo-incompressible mode, the density fluctuations are defined as the product of mass-weighted potential temperature ``P`` and the fluctuations of the inverse potential temperature ``\\chi'``. Since ``P`` is constant, this is reduced to the difference between the full density and the background density, i.e.

```math
\\rho' = P \\chi' = \\frac{P}{\\theta} - \\frac{P}{\\overline{\\theta}} = \\rho - \\overline{\\rho}.
```

# Arguments

  - `state`: Model state.
  - `model`: Dynamic equations.
"""
function synchronize_density_fluctuations!(
    state::State,
    model::PseudoIncompressible,
)
    (; rho, rhop) = state.variables.predictands

    rhop .= rho

    return
end

"""
```julia
synchronize_density_fluctuations!(state::State, model::Compressible)
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho`.

In compressible mode, the density fluctuations are defined as the product of mass-weighted potential temperature ``P`` and the fluctuations of the inverse potential temperature ``\\chi'``. In contrast to pseudo-incompressible mode, ``P`` is time-dependent, so that this is not reduced to the difference between the full density and the background density, i.e.

```math
\\rho' = P \\chi' = \\frac{P}{\\theta} - \\frac{P}{\\overline{\\theta}} = \\rho - \\frac{P}{\\overline{\\theta}}.
```

# Arguments

  - `state`: Model state.
  - `model`: Dynamic equations.
"""
function synchronize_density_fluctuations!(state::State, model::Compressible)
    (; rhostrattfc, thetastrattfc, pstrattfc) = state.atmosphere
    (; rho, rhop) = state.variables.predictands

    rhop .= rho .+ rhostrattfc .- pstrattfc ./ thetastrattfc

    return
end
