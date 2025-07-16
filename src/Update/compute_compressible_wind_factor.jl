"""
```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
```

Model-dispatched wind factor computation for momentum equation scaling.

Applied in momentum updates as: `u_new = u_old + dt * force * wind_factor`
where the wind factor ensures proper compressible flow dynamics.

# Arguments

  - `state`: Simulation state
  - `indices`: Grid indices (ix, jy, kz)
  - `variable`: Variable type for which to compute factor

# Returns

  - `::AbstractFloat`: Wind factor.
"""
function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    return compute_compressible_wind_factor(state, indices, variable, model)
end

"""
```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    model::AbstractModel,
)
```

Unity wind factor for non-compressible models.

Non-compressible models use constant density assumptions, eliminating
the need for pressure-weighted momentum scaling.

# Arguments

  - `state`: Simulation state
  - `indices`: Grid indices (ix, jy, kz)
  - `variable`: Variable type for which to compute factor
  - `model`: Model type for dispatch

# Returns

  - `::AbstractFloat`: Wind factor (`1.0`).
"""
function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    model::AbstractModel,
)
    return 1.0
end

"""
```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    model::Compressible,
)
```

Pressure-weighted wind factor for zonal velocity.

Computes Jacobian-weighted pressure average between horizontally adjacent cells
to account for variable density effects on zonal momentum transport.

Averages pressure between grid points (i, i+1) weighted by grid Jacobians:
factor = (J[i]*P[i] + J[i+1]*P[i+1]) / 2

# Arguments

  - `state`: Simulation state
  - `indices`: Grid indices (ix, jy, kz)
  - `variable`: Variable type for which to compute factor
  - `model`: Model type for dispatch

# Returns

  - `::AbstractFloat`: Zonal wind factor (``\\left(J P\\right)_{i + 1 / 2}``).
"""
function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz] * p[ix, jy, kz] +
        jac[ix + 1, jy, kz] * p[ix + 1, jy, kz]
    ) / 2
end

"""
```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    model::Compressible,
)
```

Pressure-weighted wind factor for meridional velocity.

Computes Jacobian-weighted pressure average between meridionally adjacent cells
for proper meridional momentum equation scaling in compressible flows.

Averages pressure between grid points (j, j+1) weighted by grid Jacobians:
factor = (J[i]*P[i] + J[i+1]*P[i+1]) / 2

# Arguments

  - `state`: Simulation state
  - `indices`: Grid indices (ix, jy, kz)
  - `variable`: Variable type for which to compute factor
  - `model`: Model type for dispatch

# Returns

  - `::AbstractFloat`: Meridional wind factor (`\\left(J P\\right)_{j + 1 / 2}``).
"""
function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz] * p[ix, jy, kz] +
        jac[ix, jy + 1, kz] * p[ix, jy + 1, kz]
    ) / 2
end

"""
```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
```

Pressure-weighted wind factor for vertical velocity in compressible flows.

Computes double Jacobian-weighted pressure average between vertically adjacent cells
accounting for terrain-following coordinate transformations and variable density effects
on vertical momentum transport.

Formula: `J[k]*J[k+1]*(P[k] + P[k+1]) / (J[k] + J[k+1])`

# Arguments

  - `state`: Simulation state
  - `indices`: Grid indices (ix, jy, kz)
  - `variable`: Variable type for which to compute factor
  - `model`: Model type for dispatch

# Returns

  - `::AbstractFloat`: Vertical wind factor (``\\left(J P\\right)_{k + 1 / 2}``).
"""
function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return jac[ix, jy, kz] *
           jac[ix, jy, kz + 1] *
           (p[ix, jy, kz] + p[ix, jy, kz + 1]) /
           (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
end
