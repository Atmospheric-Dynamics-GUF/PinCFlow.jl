"""
```julia
synchronize_density_fluctuations!(state::State)
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho` by dispatching to a model-specific method.

```julia
synchronize_density_fluctuations!(state::State, model::Boussinesq)
```

Return in Boussinesq mode.

In Boussinesq mode, density fluctuations don't require synchronization,
since the density is assumed constant except in the buoyancy equation.

```julia
synchronize_density_fluctuations!(state::State, model::PseudoIncompressible)
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho`.

The density fluctuations are defined as the product of the mass-weighted potential temperature and the fluctuations of the inverse potential temperature. In pseudo-incompressible mode, ``P`` is constant, so that this is reduced to the difference between ``\\rho`` and ``\\overline{\\rho}``, i.e.

```math
\\rho' = \\frac{P}{\\theta} - \\frac{P}{\\overline{\\theta}} = \\rho - \\overline{\\rho}.
```

```julia
synchronize_density_fluctuations!(state::State, model::Compressible)
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho`.

In compressible mode, ``P`` is time-dependent, so that the density fluctuations are not reduced to the difference between ``\\rho`` and ``\\overline{\\rho}``, i.e.

```math
\\rho' = \\frac{P}{\\theta} - \\frac{P}{\\overline{\\theta}} = \\rho - \\frac{P}{\\overline{\\theta}}.
```

# Arguments

  - `state`: Model state.

  - `model`: Dynamic equations.
"""
function synchronize_density_fluctuations! end

function synchronize_density_fluctuations!(state::State)
    (; model) = state.namelists.setting
    synchronize_density_fluctuations!(state, model)
    return
end

function synchronize_density_fluctuations!(state::State, model::Boussinesq)
    return
end

function synchronize_density_fluctuations!(
    state::State,
    model::PseudoIncompressible,
)
    (; rho, rhop) = state.variables.predictands

    rhop .= rho

    return
end

function synchronize_density_fluctuations!(state::State, model::Compressible)
    (; rhostrattfc, thetastrattfc, pstrattfc) = state.atmosphere
    (; rho, rhop) = state.variables.predictands

    rhop .= rho .+ rhostrattfc .- pstrattfc ./ thetastrattfc

    return
end
