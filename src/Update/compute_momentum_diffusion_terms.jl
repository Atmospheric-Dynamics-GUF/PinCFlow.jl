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

Compute and return the diffusive zonal momentum fluxes in ``\\widehat{x}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} u\\right)}^{\\widehat{x}} = \\frac{u_{i + 1 / 2} - u_{i - 1 / 2}}{\\Delta \\widehat{x}} + G^{13} \\frac{u_{k + 1} - u_{k - 1}}{2\\Delta \\widehat{z}}.
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

Compute and return the diffusive zonal momentum fluxes in ``\\widehat{y}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} u\\right)}^{\\widehat{y}} = \\frac{u_{j + 1} - u_{j - 1}}{2\\Delta \\widehat{y}} + G^{23} \\frac{u_{k + 1} - u_{k - 1}}{2\\Delta \\widehat{z}}.
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

Compute and return the diffusive zonal momentum fluxes in ``\\widehat{z}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} u\\right)}^{\\widehat{z}} = G^{13}\\frac{u_{i + 1 / 2} - u_{i - 1 / 2}}{\\Delta \\widehat{x}} + G^{23} \\frac{u_{j + 1} - u_{j - 1}}{2 \\Delta \\widehat{y}} + G^{33} \\frac{u_{k + 1} - u_{k - 1}}{2 \\Delta \\widehat{z}}.
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

Compute and return the diffusive meridional momentum fluxes in ``\\widehat{x}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} v\\right)}^{\\widehat{x}} = \\frac{v_{i + 1} - v_{i - 1}}{2 \\Delta \\widehat{x}} + G^{13} \\frac{v_{k + 1} - v_{k - 1}}{2 \\Delta \\widehat{z}}.
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

Compute and return the diffusive meridional momentum fluxes in ``\\widehat{y}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} v\\right)}^{\\widehat{y}} = \\frac{v_{j + 1 / 2} - v_{j - 1 / 2}}{\\Delta \\widehat{y}} + G^{23} \\frac{v_{k + 1} - v_{k - 1}}{2 \\Delta \\widehat{z}}.
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

Compute and return the diffusive meridional momentum fluxes in ``\\widehat{z}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} v\\right)}^{\\widehat{z}} = G^{13}\\frac{v_{i + 1} - v_{i - 1}}{2 \\Delta \\widehat{x}} + G^{23} \\frac{v_{j + 1 / 2} - v_{j - 1 / 2}}{\\Delta \\widehat{y}} + G^{33} \\frac{v_{k + 1} - v_{k - 1}}{2 \\Delta \\widehat{z}}.
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

Compute and return the diffusive vertical momentum fluxes in ``\\widehat{x}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} w\\right)}^{\\widehat{x}} = \\frac{w_{i + 1} - w_{i - 1}}{2 \\Delta \\widehat{x}} + G^{13} \\frac{w_{k + 1 / 2} - w_{k - 1 / 2}}{\\Delta \\widehat{z}}.
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

Compute and return the diffusive vertical momentum fluxes in ``\\widehat{y}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} w\\right)}^{\\widehat{y}} = \\frac{w_{j + 1} - w_{j - 1}}{2 \\Delta \\widehat{y}} + G^{23} \\frac{w_{k + 1 / 2} - w_{k - 1 / 2}}{\\Delta \\widehat{z}}.
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

Compute and return the diffusive vertical momentum fluxes in ``\\widehat{z}``-direction, i.e.

```math
\\widehat{\\left(\\boldsymbol{\\nabla} w\\right)}^{\\widehat{z}} = G^{13} \\frac{w_{i + 1} - w_{i - 1}}{2 \\Delta \\widehat{x}} + G^{23} \\frac{w_{j + 1} - w_{j - 1}}{2 \\Delta \\widehat{y}} + G^{33} \\frac{w_{k + 1 / 2} - w_{k - 1 / 2}}{\\Delta \\widehat{z}}.
```

# Arguments

  - `state`: Model state.

  - `i`: Grid-cell index on the ``\\widehat{x}``-axis.

  - `j`: Grid-cell index on the ``\\widehat{y}``-axis.

  - `k`: Grid-cell index on the ``\\widehat{z}``-axis.

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

    @ivy uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    @ivy ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    @ivy diffux =
        (u[i, j, k] - u[i - 1, j, k]) / dx +
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

    @ivy uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k])
    @ivy ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k])
    @ivy uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    @ivy ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    @ivy diffuy =
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

    @ivy uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k])
    @ivy ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k])
    @ivy uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    @ivy ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    @ivy diffuz =
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

    @ivy vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
    @ivy vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k])
    @ivy vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    @ivy vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    @ivy diffvx =
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

    @ivy vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    @ivy vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    @ivy diffvy =
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

    @ivy vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
    @ivy vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k])
    @ivy vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    @ivy vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    @ivy diffvz =
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
            compute_vertical_wind(i + 1, j, k, state) +
            compute_vertical_wind(i + 1, j, k - 1, state)
        )
    wl =
        0.5 * (
            compute_vertical_wind(i - 1, j, k, state) +
            compute_vertical_wind(i - 1, j, k - 1, state)
        )

    @ivy diffwx =
        (wr - wl) / (2.0 * dx) +
        met[i, j, k, 1, 3] * (
            compute_vertical_wind(i, j, k, state) -
            compute_vertical_wind(i, j, k - 1, state)
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
            compute_vertical_wind(i, j + 1, k, state) +
            compute_vertical_wind(i, j + 1, k - 1, state)
        )
    wb =
        0.5 * (
            compute_vertical_wind(i, j - 1, k, state) +
            compute_vertical_wind(i, j - 1, k - 1, state)
        )

    @ivy diffwy =
        (wf - wb) / (2.0 * dy) +
        met[i, j, k, 2, 3] * (
            compute_vertical_wind(i, j, k, state) -
            compute_vertical_wind(i, j, k - 1, state)
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

    @ivy diffwz =
        met[i, j, k, 1, 3] * (wr - wl) / (2.0 * dx) +
        met[i, j, k, 2, 3] * (wf - wb) / (2.0 * dy) +
        met[i, j, k, 3, 3] * (
            compute_vertical_wind(i, j, k, state) -
            compute_vertical_wind(i, j, k - 1, state)
        ) / dz

    return diffwz
end

function compute_momentum_diffusion_terms(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Z,
)::AbstractFloat
    (; u) = p0
    (; dx, dy, dz, met) = state.grid

    @ivy uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k])
    @ivy ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k])
    @ivy uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    @ivy ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])

    @ivy diffuz =
        met[i, j, k, 1, 3] * (u[i, j, k] - u[i - 1, j, k]) / dx +
        met[i, j, k, 2, 3] * (uf - ub) / (2.0 * dy) +
        met[i, j, k, 3, 3] * (uu - ud) / (2.0 * dz)

    return diffuz
end

function compute_momentum_diffusion_terms(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Z,
)::AbstractFloat
    (; v) = p0
    (; dx, dy, dz, met) = state.grid

    @ivy vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
    @ivy vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k])
    @ivy vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    @ivy vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])

    @ivy diffvz =
        met[i, j, k, 1, 3] * (vr - vl) / (2 * dx) +
        met[i, j, k, 2, 3] * (v[i, j, k] - v[i, j - 1, k]) / dy +
        met[i, j, k, 3, 3] * (vu - vd) / (2.0 * dz)

    return diffvz
end
