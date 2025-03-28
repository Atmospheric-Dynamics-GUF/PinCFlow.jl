struct Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
  phi::A
end

function Auxiliaries(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the auxiliaries.
  phi = OffsetArray(
    zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
    (-nbx):(nx + nbx),
    (-nby):(ny + nby),
    (-nbz):(nz + nbz),
  )

  # Return an Auxiliaries instance.
  return Auxiliaries(phi)
end
