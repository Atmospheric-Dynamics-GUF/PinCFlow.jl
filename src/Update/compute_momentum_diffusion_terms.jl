"""
```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::X,
)::AbstractFloat
```

Compute and return the diffusive zonal momentum fluxes in x-direction.

```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Y,
)::AbstractFloat
```

Compute and return the diffusive zonal momentum fluxes in y-direction.

```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Z,
)::AbstractFloat
```

Compute and return the diffusive zonal momentum fluxes in z-direction.

```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::X,
)::AbstractFloat
```

Compute and return the diffusive meridional momentum fluxes in x-direction.

```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Y,
)::AbstractFloat
```

Compute and return the diffusive meridional momentum fluxes in y-direction.

```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Z,
)::AbstractFloat
```

Compute and return the diffusive meridional momentum fluxes in z-direction.

```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::X,
)::AbstractFloat
```

Compute and return the diffusive vertical momentum fluxes in x-direction.

```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::Y,
)::AbstractFloat
```

Compute and return the diffusive vertical momentum fluxes in y-direction.

```julia
compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::Z,
)::AbstractFloat
```

Compute and return the diffusive vertical momentum fluxes in z-direction.

# Arguments 

  - `state`: Model state.

  - `i`: Grid-cell index on the x-axis.

  - `j`: Grid-cell index on the y-axis.

  - `k`: Grid-cell index on the z-axis. 

  - `variable`: Wind direction. 

  - `direction`: Direction of the flux.

"""

function compute_momentum_diffusion_terms end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::X,
)::AbstractFloat
    (; u) = state.variables.predictands
    (; dx, dz, met) = state.grid

    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    diffux =
        (u[i, j, k] - u[i - i, j, k]) / dx +
        met[i, j, k, 1, 3] * (uu - ud) / (2.0 * dz)

    return diffux
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Y,
)::AbstractFloat
    (; u) = state.variables.predictands
    (; dy, dz, met) = state.grid

    uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k])
    ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k])
    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    diffuy =
        (uf - ub) / (2.0 * dy) + met[i, j, k, 2, 3] * (uu - ud) / (2.0 * dz)

    return diffuy
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Z,
)::AbstractFloat
    (; u) = state.variables.predictands
    (; dx, dy, dz, met) = state.grid

    uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k])
    ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k])
    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    diffuz =
        met[i, j, k, 1, 3] * (u[i, j, k] - u[i - 1, j, k]) / dx +
        met[i, j, k, 2, 3] * (uf - ub) / (2.0 * dy) +
        met[i, j, k, 3, 3] * (uu - ud) / (2.0 * dz)

    return diffuz
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::X,
)::AbstractFloat
    (; v) = state.variables.predictands
    (; dx, dz, met) = state.grid

    vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
    vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k])
    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    diffvx =
        (vr - vl) / (2.0 * dx) + met[i, j, k, 1, 3] * (vu - vd) / (2.0 * dz)

    return diffvx
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Y,
)::AbstractFloat
    (; v) = state.variables.predictands
    (; dy, dz, met) = state.grid

    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    diffvy =
        (v[i, j, k] - v[i, j - 1, k]) / dy +
        met[i, j, k, 2, 3] * (vu - vd) / (2.0 * dz)

    return diffvy
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Z,
)::AbstractFloat
    (; v) = state.variables.predictands
    (; dx, dy, dz, met) = state.grid

    vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
    vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k])
    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    diffvz =
        met[i, j, k, 1, 3] * (vr - vl) / (2 * dx) +
        met[i, j, k, 2, 3] * (v[i, j, k] - v[i, j - 1, k]) / dy +
        met[i, j, k, 3, 3] * (vu - vd) / (2.0 * dz)

    return diffvz
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::X,
)::AbstractFloat
    (; dx, dz, met) = state.grid
    (; predictands) = state.variables
    (; grid) = state

    wr =
        0.5 * (
            compute_vertical_wind(i + 1, j, k, predictands, grid) +
            compute_vertical_wind(i + 1, j, k - 1, predictands, grid)
        )
    wl =
        0.5 * (
            compute_vertical_wind(i - 1, j, k, predictands, grid) +
            compute_vertical_wind(i - 1, j, k - 1, predictands, grid)
        )

    diffwx =
        (wr - wl) / (2.0 * dx) +
        met[i, j, k, 1, 3] * (
            compute_vertical_wind(i, j, k, predictands, grid) -
            compute_vertical_wind(i, j, k - 1, predictands, grid)
        ) / dz

    return diffwx
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::Y,
)::AbstractFloat
    (; dy, dz, met) = state.grid
    (; predictands) = state.variables
    (; grid) = state

    wf =
        0.5 * (
            compute_vertical_wind(i, j + 1, k, predictands, grid) +
            compute_vertical_wind(i, j + 1, k - 1, predictands, grid)
        )
    wb =
        0.5 * (
            compute_vertical_wind(i, j - 1, k, predictands, grid) +
            compute_vertical_wind(i, j - 1, k - 1, predictands, grid)
        )

    diffwy =
        (wf - wb) / (2.0 * dy) +
        met[i, j, k, 2, 3] * (
            compute_vertical_wind(i, j, k, predictands, grid) -
            compute_vertical_wind(i, j, k - 1, predictands, grid)
        ) / dz

    return diffwy
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::Z,
)::AbstractFloat
    (; dx, dy, dz, met) = state.grid
    (; predictands) = state.variables
    (; grid) = state

    wr =
        0.5 * (
            compute_vertical_wind(i + 1, j, k, predictands, grid) +
            compute_vertical_wind(i + 1, j, k - 1, predictands, grid)
        )
    wl =
        0.5 * (
            compute_vertical_wind(i - 1, j, k, predictands, grid) +
            compute_vertical_wind(i - 1, j, k - 1, predictands, grid)
        )
    wf =
        0.5 * (
            compute_vertical_wind(i, j + 1, k, predictands, grid) +
            compute_vertical_wind(i, j + 1, k - 1, predictands, grid)
        )
    wb =
        0.5 * (
            compute_vertical_wind(i, j - 1, k, predictands, grid) +
            compute_vertical_wind(i, j - 1, k - 1, predictands, grid)
        )

    diffwz =
        met[i, j, k, 1, 3] * (wr - wl) / (2.0 * dx) +
        met[i, j, k, 2, 3] * (wf - wb) / (2.0 * dy) +
        met[i, j, k, 3, 3] * (
            compute_vertical_wind(i, j, k, predictands, grid) -
            compute_vertical_wind(i, j, k - 1, predictands, grid)
        ) / dz

    return diffwz
end
