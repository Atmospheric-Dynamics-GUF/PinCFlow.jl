"""
```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
```

Model-dispatched wind factor computation for momentum equation scaling.

The compressible wind factor accounts for pressure-density coupling in momentum equations.
In compressible flows, momentum equations require pressure-weighted scaling to maintain
proper balance between kinetic and potential energy components.

# Arguments

  - `state::State`: Simulation state
  - `indices::NTuple{3, <:Integer}`: Grid indices (ix, jy, kz)
  - `variable::AbstractVariable`: Variable type for which to compute factor

# Returns

  - `AbstractFloat`: Pressure-weighted scaling factor

# Usage

Applied in momentum updates as: `u_new = u_old + dt * force * wind_factor`
where the wind factor ensures proper compressible flow dynamics.
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

# Returns

  - `Float64`: Always returns 1.0
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

# Implementation

Averages pressure between grid points (i, i+1) weighted by grid Jacobians:
factor = (J[i]*P[i] + J[i+1]*P[i+1]) / 2

# Returns

  - `AbstractFloat`: Pressure-weighted factor for zonal momentum scaling
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

# Implementation

Averages pressure between grid points (j, j+1) weighted by grid Jacobians:
factor = (J[i]*P[i] + J[i+1]*P[i+1]) / 2

# Returns

  - `AbstractFloat`: Pressure-weighted factor for meridional momentum scaling
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

# Returns

  - `AbstractFloat`: Pressure-weighted scaling factor

# Implementation

Formula: `J[k]*J[k+1]*(P[k] + P[k+1]) / (J[k] + J[k+1])`
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
