abstract type AbstractCoordinate end
struct Cartesian <: AbstractCoordinate end
struct TFC <: AbstractCoordinate end

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
      jac[i, j, k + 1] * (met[i, j, k, 1, 3] * uc + met[i, j, k, 2, 3] * vc) +
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
    jac[i, j, k] * (met[i, j, k + 1, 1, 3] * uu + met[i, j, k + 1, 2, 3] * vu)
  ) / (jac[i, j, k] + jac[i, j, k + 1]) + wedgeu / jacedgeu
end
