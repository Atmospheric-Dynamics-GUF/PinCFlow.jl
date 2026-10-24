"""
```julia
synchronize_density_fluctuations!(state::State)::Nothing
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho` by dispatching to a model-specific method.

```julia
synchronize_density_fluctuations!(
    state::State,
    model::Val{:Boussinesq},
)::Nothing
```

Return in Boussinesq mode.

In Boussinesq mode, density fluctuations don't require synchronization,
since the density is assumed constant except in the buoyancy equation.

```julia
synchronize_density_fluctuations!(
    state::State,
    model::Val{:PseudoIncompressible},
)::Nothing
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho`.

The density fluctuations are defined as the product of the mass-weighted potential temperature and the fluctuations of the inverse potential temperature. In pseudo-incompressible mode, ``P`` is constant, so that this is reduced to the difference between ``\\rho`` and ``\\bar{\\rho}``, i.e.

```math
\\rho' = \\frac{P}{\\theta} - \\frac{P}{\\bar{\\theta}} = \\rho - \\bar{\\rho}.
```

```julia
synchronize_density_fluctuations!(
    state::State,
    model::Val{:Compressible},
)::Nothing
```

Synchronize the density fluctuations in `state.variables.predictands.rhop` with the density in `state.variables.predictands.rho`.

In compressible mode, ``P`` is time-dependent, so that the density fluctuations are not reduced to the difference between ``\\rho`` and ``\\bar{\\rho}``, i.e.

```math
\\rho' = \\frac{P}{\\theta} - \\frac{P}{\\bar{\\theta}} = \\rho - \\frac{P}{\\bar{\\theta}}.
```

# Arguments

  - `state`: Model state.

  - `model`: Dynamic equations.
"""
function synchronize_density_fluctuations! end

function synchronize_density_fluctuations!(state::State)::Nothing
    (; model) = state.namelists.atmosphere
    @dispatch_model synchronize_density_fluctuations!(state, Val(model))
    nothing
end

function synchronize_density_fluctuations!(
    state::State,
    model::Val{:Boussinesq},
)::Nothing
    nothing
end

function synchronize_density_fluctuations!(
    state::State,
    model::Val{:PseudoIncompressible},
)::Nothing
    (; rho, rhop) = state.variables.predictands

    rhop .= rho

    nothing
end

function synchronize_density_fluctuations!(
    state::State,
    model::Val{:Compressible},
)::Nothing
    (; rhobar, thetabar, pbar) = state.atmosphere
    (; rho, rhop) = state.variables.predictands

    rhop .= rho .+ rhobar .- pbar ./ thetabar

    nothing
end
