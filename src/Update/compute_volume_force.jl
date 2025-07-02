"""
```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
```

Test case-dispatched volume force computation.

Dispatches to specific implementation based on test case type.

# Arguments

  - `state::State`: Simulation state
  - `indices::NTuple{3, <:Integer}`: Grid indices (ix, jy, kz)
  - `variable::AbstractVariable`: Variable type for force computation

# Returns

  - `AbstractFloat`: Volume force per unit mass
"""
function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; testcase) = state.namelists.setting

    return compute_volume_force(state, indices, variable, testcase)
end

"""
```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    testcase::AbstractTestCase,
)
```

Default zero volume force for standard test cases.

Most test cases have no additional volume forcing beyond standard physics.

# Returns

  - `Float64`: Always returns 0.0
"""
function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    testcase::AbstractTestCase,
)
    return 0.0
end

"""
```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    testcase::AbstractWKBTestCase,
)
```

Zonal volume force from gravity wave drag.

Computes cell-edge interpolated zonal wind tendency from WKB gravity wave parameterization for use in momentum equations.

# Implementation

Averages `dudt` between horizontally adjacent cells:
`force = (dudt[i] + dudt[i+1]) / 2`

# Returns

  - `AbstractFloat`: Zonal acceleration [m/s²]
"""
function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    testcase::AbstractWKBTestCase,
)
    (; dudt) = state.wkb.tendencies
    (ix, jy, kz) = indices

    return (dudt[ix, jy, kz] + dudt[ix + 1, jy, kz]) / 2
end

"""
```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    testcase::AbstractWKBTestCase,
)
```

Meridional volume force from gravity wave drag.

Computes cell-edge interpolated meridional wind tendency for momentum equations.

# Implementation

Averages `dvdt` between meridionally adjacent cells.

# Returns

  - `AbstractFloat`: Meridional acceleration [m/s²]
"""
function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    testcase::AbstractWKBTestCase,
)
    (; dvdt) = state.wkb.tendencies
    (ix, jy, kz) = indices

    return (dvdt[ix, jy, kz] + dvdt[ix, jy + 1, kz]) / 2
end

"""
```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    testcase::AbstractWKBTestCase,
)
```

Vertical volume force from terrain-following coordinate transformation.

Transforms horizontal gravity wave tendencies to vertical component using
metric tensor coefficients and Jacobian weighting for terrain-following coordinates.

# Implementation

  - **Metric transformation**: Uses `met` tensor components
  - **Jacobian weighting**: Properly averages between vertical levels
  - **Coordinate mapping**: `w_force = (∂ζ/∂x)·u_force + (∂ζ/∂y)·v_force`

# Returns

  - `AbstractFloat`: Vertical acceleration [m/s²]
"""
function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    testcase::AbstractWKBTestCase,
)
    (; jac, met) = state.grid
    (; dudt, dvdt) = state.wkb.tendencies
    (ix, jy, kz) = indices

    return (
        jac[ix, jy, kz + 1] * (
            met[ix, jy, kz, 1, 3] * dudt[ix, jy, kz] +
            met[ix, jy, kz, 2, 3] * dvdt[ix, jy, kz]
        ) +
        jac[ix, jy, kz] * (
            met[ix, jy, kz + 1, 1, 3] * dudt[ix, jy, kz + 1] +
            met[ix, jy, kz + 1, 2, 3] * dvdt[ix, jy, kz + 1]
        )
    ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
end

"""
```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::P,
    testcase::AbstractWKBTestCase,
)
```

Thermal forcing from gravity wave heating.

Returns potential temperature tendency from gravity wave energy dissipation
for use in thermodynamic equation.

# Returns

  - `AbstractFloat`: Heating rate [K/s]
"""
function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::P,
    testcase::AbstractWKBTestCase,
)
    (; dthetadt) = state.wkb.tendencies
    (ix, jy, kz) = indices

    return dthetadt[ix, jy, kz]
end
