"""
```julia
transform(
    i::Integer,
    j::Integer,
    k::Integer,
    uedger::AbstractFloat,
    uuedger::AbstractFloat,
    uedgel::AbstractFloat,
    uuedgel::AbstractFloat,
    vedgef::AbstractFloat,
    vuedgef::AbstractFloat,
    vedgeb::AbstractFloat,
    vuedgeb::AbstractFloat,
    wedgeu::AbstractFloat,
    coordinate::Cartesian,
    state::State,
)::AbstractFloat
```

Perform the transformation of a vertical-wind-like variable from the transformed system to the Cartesian one, given the wind-like components at the grid points surrounding ``\\left(i, j, k + 1 / 2\\right)``, and return the result.

The discretized transformation rule for the vertical wind is given by

```math
w_{k + 1 / 2} = J_{k + 1 / 2} \\left[- \\left(G^{1 3} u\\right)_{k + 1 / 2} - \\left(G^{2 3} v\\right)_{k + 1 / 2} + \\widehat{w}_{k + 1 / 2}\\right].
```

```julia
transform(
    i::Integer,
    j::Integer,
    k::Integer,
    uedger::AbstractFloat,
    uuedger::AbstractFloat,
    uedgel::AbstractFloat,
    uuedgel::AbstractFloat,
    vedgef::AbstractFloat,
    vuedgef::AbstractFloat,
    vedgeb::AbstractFloat,
    vuedgeb::AbstractFloat,
    wedgeu::AbstractFloat,
    coordinate::Transformed,
    state::State,
)::AbstractFloat
```

Perform the transformation of a vertical-wind-like variable from the Cartesian system to the transformed one, given the wind-like components at the grid points surrounding ``\\left(i, j, k + 1 / 2\\right)``, and return the result.

The discretized transformation rule for the vertical wind is given by

```math
\\widehat{w}_{k + 1 / 2} = \\left(G^{1 3} u\\right)_{k + 1 / 2} + \\left(G^{2 3} v\\right)_{k + 1 / 2} + \\frac{w_{k + 1 / 2}}{J_{k + 1 / 2}}.
```

# Arguments

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `uedger`: Zonal-wind equivalent at ``\\left(i + 1 / 2, j, k\\right)``.

  - `uuedger`: Zonal-wind equivalent at ``\\left(i + 1 / 2, j, k + 1\\right)``.

  - `uedgel`: Zonal-wind equivalent at ``\\left(i - 1 / 2, j, k\\right)``.

  - `uuedgel`: Zonal-wind equivalent at ``\\left(i - 1 / 2, j, k + 1\\right)``.

  - `vedgef`: Meridional-wind equivalent at ``\\left(i, j + 1 / 2, k\\right)``.

  - `vuedgef`: Meridional-wind equivalent at ``\\left(i, j + 1 / 2, k + 1\\right)``.

  - `vedgeb`: Meridional-wind equivalent at ``\\left(i, j - 1 / 2, k\\right)``.

  - `vuedgeb`: Meridional-wind equivalent at ``\\left(i, j - 1 / 2, k + 1\\right)``.

  - `wedgeu`: Transformed-vertical-wind equivalent at ``\\left(i, j, k + 1 / 2\\right)``

  - `coordinate`: Coordinate system to transform to.

  - `state`: Model state.
"""
function transform end

function transform(
    i::Integer,
    j::Integer,
    k::Integer,
    uedger::AbstractFloat,
    uuedger::AbstractFloat,
    uedgel::AbstractFloat,
    uuedgel::AbstractFloat,
    vedgef::AbstractFloat,
    vuedgef::AbstractFloat,
    vedgeb::AbstractFloat,
    vuedgeb::AbstractFloat,
    wedgeu::AbstractFloat,
    coordinate::Cartesian,
    state::State,
)::AbstractFloat
    (; jac, met) = state.grid

    @ivy jacedgeu =
        2.0 * jac[i, j, k] * jac[i, j, k + 1] /
        (jac[i, j, k] + jac[i, j, k + 1])

    @ivy uc = 0.5 * (uedger + uedgel)
    @ivy uu = 0.5 * (uuedger + uuedgel)
    @ivy vc = 0.5 * (vedgef + vedgeb)
    @ivy vu = 0.5 * (vuedgef + vuedgeb)

    @ivy return jacedgeu * (
        -(
            jac[i, j, k + 1] *
            (met[i, j, k, 1, 3] * uc + met[i, j, k, 2, 3] * vc) +
            jac[i, j, k] *
            (met[i, j, k + 1, 1, 3] * uu + met[i, j, k + 1, 2, 3] * vu)
        ) / (jac[i, j, k] + jac[i, j, k + 1]) + wedgeu
    )
end

function transform(
    i::Integer,
    j::Integer,
    k::Integer,
    uedger::AbstractFloat,
    uuedger::AbstractFloat,
    uedgel::AbstractFloat,
    uuedgel::AbstractFloat,
    vedgef::AbstractFloat,
    vuedgef::AbstractFloat,
    vedgeb::AbstractFloat,
    vuedgeb::AbstractFloat,
    wedgeu::AbstractFloat,
    coordinate::Transformed,
    state::State,
)::AbstractFloat
    (; jac, met) = state.grid

    @ivy jacedgeu =
        2.0 * jac[i, j, k] * jac[i, j, k + 1] /
        (jac[i, j, k] + jac[i, j, k + 1])

    @ivy uc = 0.5 * (uedger + uedgel)
    @ivy uu = 0.5 * (uuedger + uuedgel)
    @ivy vc = 0.5 * (vedgef + vedgeb)
    @ivy vu = 0.5 * (vuedgef + vuedgeb)

    @ivy return (
        jac[i, j, k + 1] * (met[i, j, k, 1, 3] * uc + met[i, j, k, 2, 3] * vc) +
        jac[i, j, k] *
        (met[i, j, k + 1, 1, 3] * uu + met[i, j, k + 1, 2, 3] * vu)
    ) / (jac[i, j, k] + jac[i, j, k + 1]) + wedgeu / jacedgeu
end
