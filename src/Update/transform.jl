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
)
```

Transform vertical velocity from Cartesian to terrain-following coordinates.

# Coordinate Transformation

  - **Cartesian**: `w^ζ = w - (g^{1ζ}u + g^{2ζ}v)`
  - **Metric terms**: Uses `met[i,j,k,1,3]` and `met[i,j,k,2,3]`
  - **Jacobian weighting**: Harmonic mean between adjacent cells

# Arguments

  - `i, j, k`: Grid indices
  - `u*, v*, w*`: Velocity components at cell edges
  - `coordinate::Cartesian`: Coordinate system type
  - `grid`: Structure containing metric tensor and Jacobian

# Returns

  - `AbstractFloat`: Contravariant vertical velocity w^ζ
"""
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
)
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
    coordinate::Transformed,
    grid::Grid,
)
```

Transform vertical velocity including terrain-following coordinate effects.

For terrain-following coordinates with additional metric corrections.

# Arguments

  - `i, j, k`: Grid indices
  - `u*, v*, w*`: Velocity components at cell edges
  - `coordinate::Transformed`: Type dispatch for transformed coordinate system
  - `grid`: Structure containing metric tensor and Jacobian

# Returns

  - `AbstractFloat`: Contravariant vertical velocity w^ζ
"""
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
)
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
