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
    2.0 * jac[i, j, k] * jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1])

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
  coordinate::TFC,
  grid::Grid,
)
  (; jac, met) = grid

  jacedgeu =
    2.0 * jac[i, j, k] * jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1])

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

function compute_vertical_wind(
  i::Integer,
  j::Integer,
  k::Integer,
  predictands::Predictands,
  grid::Grid,
)
  (; u, v, w) = predictands

  uedger = u[i, j, k]
  uuedger = u[i, j, k + 1]
  uedgel = u[i - 1, j, k]
  uuedgel = u[i - 1, j, k + 1]
  vedgef = v[i, j, k]
  vuedgef = v[i, j, k + 1]
  vedgeb = v[i, j - 1, k]
  vuedgeb = v[i, j - 1, k + 1]
  wedgeu = w[i, j, k]

  return transform(
    i,
    j,
    k,
    uedger,
    uuedger,
    uedgel,
    uuedgel,
    vedgef,
    vuedgef,
    vedgeb,
    vuedgeb,
    wedgeu,
    Cartesian(),
    grid,
  )
end