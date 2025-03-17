struct Correction{A <: OffsetArray{<:AbstractFloat, 3}}
  corx::A
  cory::B
end

function Correction(namelists::Namelists, domain::Domain)

  # Get all necessary fields.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize correction fields.
  (corx, cory) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:2
  )

  # Return a Correction instance.
  return Correction(corx, cory)
end
