struct Operator{A <: OffsetArray{<:AbstractFloat, 3}}
  s::A
end

function Operator(namelists::Namelists, domain::Domain)

  # Get all necessary fields.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize s.
  s = OffsetArray(
    zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
    (-nbx):(nx + nbx),
    (-nby):(ny + nby),
    (-nbz):(nz + nbz),
  )

  # Return an Operator instance.
  return Operator(s)
end
