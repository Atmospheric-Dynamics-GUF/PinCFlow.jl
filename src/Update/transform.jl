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
    grid::Grid,
)::AbstractFloat
```

Perform the transformation of a vertical-wind-like variable from the transformed system to the Cartesian one, given the wind-like components at the grid points surrounding `(i, j, k + 1 / 2)`, and return the result.

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
    grid::Grid,
)::AbstractFloat
```

Perform the transformation of a vertical-wind-like variable from the Cartesian system to the transformed one, given the wind-like components at the grid points surrounding `(i, j, k + 1 / 2)`, and return the result.

The discretized transformation rule for the vertical wind is given by

```math
\\widehat{w}_{k + 1 / 2} = \\left(G^{1 3} u\\right)_{k + 1 / 2} + \\left(G^{2 3} v\\right)_{k + 1 / 2} + \\frac{w_{k + 1 / 2}}{J_{k + 1 / 2}}.
```

# Arguments

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `uedger`: Zonal-wind equivalent at `(i + 1 / 2, j, k)`.

  - `uuedger`: Zonal-wind equivalent at `(i + 1 / 2, j, k + 1)`.

  - `uedgel`: Zonal-wind equivalent at `(i - 1 / 2, j, k)`.

  - `uuedgel`: Zonal-wind equivalent at `(i - 1 / 2, j, k + 1)`.

  - `vedgef`: Meridional-wind equivalent at `(i, j + 1 / 2, k)`.

  - `vuedgef`: Meridional-wind equivalent at `(i, j + 1 / 2, k + 1)`.

  - `vedgeb`: Meridional-wind equivalent at `(i, j - 1 / 2, k)`.

  - `vuedgeb`: Meridional-wind equivalent at `(i, j - 1 / 2, k + 1)`.

  - `wedgeu`: Transformed-vertical-wind equivalent at `(i, j, k + 1 / 2)`

  - `coordinate`: Coordinate system to transform to.

  - `grid`: Collection of parameters and fields that describe the grid.
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
    grid::Grid,
)::AbstractFloat
    (; jac, met) = grid

    jacedgeu =
        2.0 * jac[i, j, k] * jac[i, j, k + 1] /
        (jac[i, j, k] + jac[i, j, k + 1])

    uc = 0.5 * (uedger + uedgel)
    uu = 0.5 * (uuedger + uuedgel)
    vc = 0.5 * (vedgef + vedgeb)
    vu = 0.5 * (vuedgef + vuedgeb)

    return jacedgeu * (
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
    grid::Grid,
)::AbstractFloat
    (; jac, met) = grid

    jacedgeu =
        2.0 * jac[i, j, k] * jac[i, j, k + 1] /
        (jac[i, j, k] + jac[i, j, k + 1])

    uc = 0.5 * (uedger + uedgel)
    uu = 0.5 * (uuedger + uuedgel)
    vc = 0.5 * (vedgef + vedgeb)
    vu = 0.5 * (vuedgef + vuedgeb)

    return (
        jac[i, j, k + 1] * (met[i, j, k, 1, 3] * uc + met[i, j, k, 2, 3] * vc) +
        jac[i, j, k] *
        (met[i, j, k + 1, 1, 3] * uu + met[i, j, k + 1, 2, 3] * vu)
    ) / (jac[i, j, k] + jac[i, j, k + 1]) + wedgeu / jacedgeu
end
