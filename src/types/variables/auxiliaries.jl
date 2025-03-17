struct Auxiliaries{A <: OffsetArray{<:AbstractFloat, 3}}
  rhobar::A
  rhopbar::A
  ubar::A
  vbar::A
  wbar::A
end

function Auxiliaries(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the auxiliaries.
  (rhobar, rhopbar, ubar, vbar, wbar) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:5
  )

  # Return an Auxiliaries instance.
  return Auxiliaries(rhobar, rhopbar, ubar, vbar, wbar)
end
