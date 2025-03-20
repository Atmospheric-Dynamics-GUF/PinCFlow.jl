struct Backups{A <: AbstractArray{<:AbstractFloat, 3}}
  rhoold::A
  rhopold::A
  uold::A
  vold::A
  wold::A
end

function Backups(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the backups.
  (rhoold, rhopold, uold, vold, wold) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:5
  )

  # Return a Backups instance.
  return Backups(rhoold, rhopold, uold, vold, wold)
end
