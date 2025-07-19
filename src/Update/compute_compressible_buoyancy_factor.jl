"""
```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
```

Compute the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` or ``\\left(i, j, k + 1 / 2\\right)``, depending on `variable`.

In pseudo-incompressible mode, the squared buoyancy frequency used by PinCFlow is

```math
N^2 = \\frac{g}{\\overline{\\theta}} \\frac{\\mathrm{d} \\overline{\\theta}}{\\mathrm{d} z},
```

whereas in compressible mode, it is

```math
N^2 = \\frac{g P}{\\rho \\overline{\\theta}^2} \\frac{\\mathrm{d} \\overline{\\theta}}{\\mathrm{d} z}.
```

In both modes, the buoyancy term is expressed in terms of ``N^2``. Thus, one has

```math
\\left(\\frac{\\partial b'}{\\partial t}\\right)_{N^2} = f_{b'} N^2 w,
```

with ``f_{b'} = 1`` in compressible mode and ``f_{b'} = \\overline{\\rho} / \\rho`` in pseudo-incompressible mode. This method returns ``f_{b'}``, either at ``\\left(i, j, k\\right)`` or at ``\\left(i, j, k + 1 / 2\\right)``, based on the type of `variable`. Note that both values of ``f_{b'}`` are equivalent in Boussinesq mode, where ``\\rho = \\overline{\\rho} = \\rho_0`` (since the density fluctuations are treated separately).

# Arguments

  - `state`: Model state.
  - `indices`: Grid-cell indices.
  - `variable`: Variable for which the factor is needed.

# Returns

  - `::AbstractFloat`: Buoyancy factor.
"""
function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    return compute_compressible_buoyancy_factor(state, indices, variable, model)
end

"""
```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::Compressible,
)
```

Return `1.0` as the factor by which the buoyancy term should be multiplied in compressible mode.

In compressible mode, the squared buoyancy frequency used by PinCFlow is

```math
N^2 = \\frac{g P}{\\rho \\overline{\\theta}^2} \\frac{\\mathrm{d} \\overline{\\theta}}{\\mathrm{d} z},
```

so that the buoyancy term may be written as

```math
\\left(\\frac{\\partial b'}{\\partial t}\\right)_{N^2} = f_{b'} N^2 w,
```

with ``f_{b'} = 1``.

# Arguments

  - `state`: Model state.
  - `indices`: Grid-cell indices.
  - `variable`: Variable for which the factor is needed.
  - `model`: Dynamic equations.

# Returns

  - `::AbstractFloat`: Buoyancy factor (`1.0`).
"""
function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::Compressible,
)
    return 1.0
end

"""
```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::AbstractModel,
)
```

Compute the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` in pseudo-incompressible mode (this method is also used in Boussinesq mode).

In pseudo-incompressible mode, the squared buoyancy frequency used by PinCFlow is

```math
N^2 = \\frac{g}{\\overline{\\theta}} \\frac{\\mathrm{d} \\overline{\\theta}}{\\mathrm{d} z},
```

so that the buoyancy term may be written as

```math
\\left(\\frac{\\partial b'}{\\partial t}\\right)_{N^2} = f_{b'} N^2 w,
```

with ``f_{b'} = \\overline{\\rho} / \\rho``. In Boussinesq mode, this is reduced to ``f_{b'} = 1`` (since ``\\rho = \\overline{\\rho} = \\rho_0``, with the density fluctuations being treated separately).

# Arguments

  - `state`: Model state.
  - `indices`: Grid-cell indices.
  - `variable`: Variable for which the factor is needed.
  - `model`: Dynamic equations.

# Returns

  - `::AbstractFloat`: Buoyancy factor (``\\overline{\\rho} / \\rho``).
"""
function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::AbstractModel,
)
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return rhostrattfc[ix, jy, kz] / (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz])
end

"""
```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
```

Return `1.0` as the factor by which the buoyancy term should be multiplied in compressible mode.

In compressible mode, the squared buoyancy frequency used by PinCFlow is

```math
N^2 = \\frac{g P}{\\rho \\overline{\\theta}^2} \\frac{\\mathrm{d} \\overline{\\theta}}{\\mathrm{d} z},
```

so that the buoyancy term may be written as

```math
\\left(\\frac{\\partial b'}{\\partial t}\\right)_{N^2} = f_{b'} N^2 w,
```

with ``f_{b'} = 1``.

# Arguments

  - `state`: Model state.
  - `indices`: Grid-cell indices.
  - `variable`: Variable for which the factor is needed.
  - `model`: Dynamic equations.

# Returns

  - `::AbstractFloat`: Buoyancy factor (`1.0`).
"""
function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
    return 1.0
end

"""
```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::AbstractModel,
)
```

Compute the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k + 1 / 2\\right)`` in pseudo-incompressible mode (this method is also used in Boussinesq mode).

In pseudo-incompressible mode, the squared buoyancy frequency used by PinCFlow is

```math
N^2 = \\frac{g}{\\overline{\\theta}} \\frac{\\mathrm{d} \\overline{\\theta}}{\\mathrm{d} z},
```

so that the buoyancy term may be written as

```math
\\left(\\frac{\\partial b'}{\\partial t}\\right)_{N^2} = f_{b'} N^2 w,
```

with ``f_{b'} = \\overline{\\rho} / \\rho``. In Boussinesq mode, this is reduced to ``f_{b'} = 1`` (since ``\\rho = \\overline{\\rho} = \\rho_0``, with the density fluctuations being treated separately). The interpolation to ``\\left(i, j, k + 1 / 2\\right)`` is

```math
f_{w'} = \\frac{\\overline{\\rho}_{k + 1 / 2}}{\\rho_{k + 1 / 2}} = \\frac{J_{k + 1} \\overline{\\rho} + J \\overline{\\rho}_{k + 1}}{J_{k + 1} \\rho + J \\rho_{k + 1}}.
```

# Arguments

  - `state`: Model state.
  - `indices`: Grid-cell indices.
  - `variable`: Variable for which the factor is needed.
  - `model`: Dynamic equations.

# Returns

  - `::AbstractFloat`: Buoyancy factor (``\\overline{\\rho}_{k + 1 / 2} / \\rho_{k + 1 / 2}``).
"""
function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::AbstractModel,
)
    (; jac) = state.grid
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz + 1] * rhostrattfc[ix, jy, kz] +
        jac[ix, jy, kz] * rhostrattfc[ix, jy, kz + 1]
    ) / (
        jac[ix, jy, kz + 1] * (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]) +
        jac[ix, jy, kz] * (rho[ix, jy, kz + 1] + rhostrattfc[ix, jy, kz + 1])
    )
end
