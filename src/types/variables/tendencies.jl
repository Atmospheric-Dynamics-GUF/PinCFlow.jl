struct Tendencies{A <: AbstractArray{<:AbstractFloat, 3}}
  drho::A
  drhop::A
  du::A
  dv::A
  dw::A
  dpip::A
end

function Tendencies(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the tendencies.
  (drho, drhop, du, dv, dw) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:5
  )

  # Initialize the pressure differences.
  dpip = OffsetArray(
    zeros((nx + 2, ny + 2, nz + 2)),
    0:(nx + 1),
    0:(ny + 1),
    0:(nz + 1),
  )

  # Return a Variables instance.
  return Tendencies(drho, drhop, du, dv, dw, dpip)
end
