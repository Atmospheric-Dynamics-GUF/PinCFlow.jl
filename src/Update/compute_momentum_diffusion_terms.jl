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

Compute and return the diffusive zonal momentum fluxes in ``\\hat{x}``-direction, i.e.

```math
\\Xi_u^{\\hat{x}} = \\frac{u_{i + 1 / 2} - u_{i - 1 / 2}}{\\Delta \\hat{x}} + G^{13} \\frac{u_{k + 1} - u_{k - 1}}{2\\Delta \\hat{z}}.
```

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

Compute and return the diffusive zonal momentum fluxes in ``\\hat{y}``-direction, i.e.

```math
\\Xi_u^{\\hat{y}} = \\frac{u_{j + 1} - u_{j - 1}}{2\\Delta \\hat{y}} + G^{23} \\frac{u_{k + 1} - u_{k - 1}}{2\\Delta \\hat{z}}.
```

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

Compute and return the diffusive zonal momentum fluxes in ``\\hat{z}``-direction, i.e.

```math
\\Xi_u^{\\hat{z}} = G^{13}\\frac{u_{i + 1 / 2} - u_{i - 1 / 2}}{\\Delta \\hat{x}} + G^{23} \\frac{u_{j + 1} - u_{j - 1}}{2 \\Delta \\hat{y}} + G^{33} \\frac{u_{k + 1} - u_{k - 1}}{2 \\Delta \\hat{z}}.
```

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

Compute and return the diffusive meridional momentum fluxes in ``\\hat{x}``-direction, i.e.

```math
\\Xi_v^{\\hat{x}} = \\frac{v_{i + 1} - v_{i - 1}}{2 \\Delta \\hat{x}} + G^{13} \\frac{v_{k + 1} - v_{k - 1}}{2 \\Delta \\hat{z}}.
```

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

Compute and return the diffusive meridional momentum fluxes in ``\\hat{y}``-direction, i.e.

```math
\\Xi_v^{\\hat{y}} = \\frac{v_{j + 1 / 2} - v_{j - 1 / 2}}{\\Delta \\hat{y}} + G^{23} \\frac{v_{k + 1} - v_{k - 1}}{2 \\Delta \\hat{z}}.
```

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

Compute and return the diffusive meridional momentum fluxes in ``\\hat{z}``-direction, i.e.

```math
\\Xi_v^{\\hat{z}} = G^{13}\\frac{v_{i + 1} - v_{i - 1}}{2 \\Delta \\hat{x}} + G^{23} \\frac{v_{j + 1 / 2} - v_{j - 1 / 2}}{\\Delta \\hat{y}} + G^{33} \\frac{v_{k + 1} - v_{k - 1}}{2 \\Delta \\hat{z}}.
```

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

Compute and return the diffusive vertical momentum fluxes in ``\\hat{x}``-direction, i.e.

```math
\\Xi_w^{\\hat{x}} = \\frac{w_{i + 1} - w_{i - 1}}{2 \\Delta \\hat{x}} + G^{13} \\frac{w_{k + 1 / 2} - w_{k - 1 / 2}}{\\Delta \\hat{z}}.
```

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

Compute and return the diffusive vertical momentum fluxes in ``\\hat{y}``-direction, i.e.

```math
\\Xi_w^{\\hat{y}} = \\frac{w_{j + 1} - w_{j - 1}}{2 \\Delta \\hat{y}} + G^{23} \\frac{w_{k + 1 / 2} - w_{k - 1 / 2}}{\\Delta \\hat{z}}.
```

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

Compute and return the diffusive vertical momentum fluxes in ``\\hat{z}``-direction, i.e.

```math
\\Xi_w^{\\hat{z}} = G^{13} \\frac{w_{i + 1} - w_{i - 1}}{2 \\Delta \\hat{x}} + G^{23} \\frac{w_{j + 1} - w_{j - 1}}{2 \\Delta \\hat{y}} + G^{33} \\frac{w_{k + 1 / 2} - w_{k - 1 / 2}}{\\Delta \\hat{z}}.
```

# Arguments

  - `state`: Model state.

  - `i`: Grid-cell index on the ``\\hat{x}``-axis.

  - `j`: Grid-cell index on the ``\\hat{y}``-axis.

  - `k`: Grid-cell index on the ``\\hat{z}``-axis.

  - `variable`: Wind direction.

  - `direction`: Direction of the flux.
"""
function compute_momentum_diffusion_terms end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::X,
)::AbstractFloat
    (; u) = state.variables.predictands
    (; dx, dz, met13, met23, met33) = state.grid

    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    diffux =
        (u[i, j, k] - u[i - 1, j, k]) / dx +
        met13[i, j, k] * (uu - ud) / (2.0 * dz)

    return diffux
end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Y,
)::AbstractFloat
    (; u) = state.variables.predictands
    (; dy, dz, met13, met23, met33) = state.grid

    uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k])
    ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k])
    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    diffuy = (uf - ub) / (2.0 * dy) + met23[i, j, k] * (uu - ud) / (2.0 * dz)

    return diffuy
end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Z,
)::AbstractFloat
    (; u) = state.variables.predictands
    (; dx, dy, dz, met13, met23, met33) = state.grid

    uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k])
    ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k])
    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    diffuz =
        met13[i, j, k] * (u[i, j, k] - u[i - 1, j, k]) / dx +
        met23[i, j, k] * (uf - ub) / (2.0 * dy) +
        met33[i, j, k] * (uu - ud) / (2.0 * dz)

    return diffuz
end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::X,
)::AbstractFloat
    (; v) = state.variables.predictands
    (; dx, dz, met13, met23, met33) = state.grid

    vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
    vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k])
    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    diffvx = (vr - vl) / (2.0 * dx) + met13[i, j, k] * (vu - vd) / (2.0 * dz)

    return diffvx
end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Y,
)::AbstractFloat
    (; v) = state.variables.predictands
    (; dy, dz, met13, met23, met33) = state.grid

    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    diffvy =
        (v[i, j, k] - v[i, j - 1, k]) / dy +
        met23[i, j, k] * (vu - vd) / (2.0 * dz)

    return diffvy
end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Z,
)::AbstractFloat
    (; v) = state.variables.predictands
    (; dx, dy, dz, met13, met23, met33) = state.grid

    vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
    vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k])
    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    diffvz =
        met13[i, j, k] * (vr - vl) / (2 * dx) +
        met23[i, j, k] * (v[i, j, k] - v[i, j - 1, k]) / dy +
        met33[i, j, k] * (vu - vd) / (2.0 * dz)

    return diffvz
end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::X,
)::AbstractFloat
    (; dx, dz, met13, met23, met33) = state.grid

    wr =
        0.5 * (
            compute_vertical_wind(i + 1, j, k, state) +
            compute_vertical_wind(i + 1, j, k - 1, state)
        )
    wl =
        0.5 * (
            compute_vertical_wind(i - 1, j, k, state) +
            compute_vertical_wind(i - 1, j, k - 1, state)
        )

    diffwx =
        (wr - wl) / (2.0 * dx) +
        met13[i, j, k] * (
            compute_vertical_wind(i, j, k, state) -
            compute_vertical_wind(i, j, k - 1, state)
        ) / dz

    return diffwx
end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::Y,
)::AbstractFloat
    (; dy, dz, met13, met23, met33) = state.grid

    wf =
        0.5 * (
            compute_vertical_wind(i, j + 1, k, state) +
            compute_vertical_wind(i, j + 1, k - 1, state)
        )
    wb =
        0.5 * (
            compute_vertical_wind(i, j - 1, k, state) +
            compute_vertical_wind(i, j - 1, k - 1, state)
        )

    diffwy =
        (wf - wb) / (2.0 * dy) +
        met23[i, j, k] * (
            compute_vertical_wind(i, j, k, state) -
            compute_vertical_wind(i, j, k - 1, state)
        ) / dz

    return diffwy
end

@ivy function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::Z,
)::AbstractFloat
    (; dx, dy, dz, met13, met23, met33) = state.grid

    wr =
        0.5 * (
            compute_vertical_wind(i + 1, j, k, state) +
            compute_vertical_wind(i + 1, j, k - 1, state)
        )
    wl =
        0.5 * (
            compute_vertical_wind(i - 1, j, k, state) +
            compute_vertical_wind(i - 1, j, k - 1, state)
        )
    wf =
        0.5 * (
            compute_vertical_wind(i, j + 1, k, state) +
            compute_vertical_wind(i, j + 1, k - 1, state)
        )
    wb =
        0.5 * (
            compute_vertical_wind(i, j - 1, k, state) +
            compute_vertical_wind(i, j - 1, k - 1, state)
        )

    diffwz =
        met13[i, j, k] * (wr - wl) / (2.0 * dx) +
        met23[i, j, k] * (wf - wb) / (2.0 * dy) +
        met33[i, j, k] * (
            compute_vertical_wind(i, j, k, state) -
            compute_vertical_wind(i, j, k - 1, state)
        ) / dz

    return diffwz
end
