struct Fluxes{A <: AbstractArray{<:AbstractFloat, 4}}
  phirho::A
  phirhop::A
  phiu::A
  phiv::A
  phiw::A
end

function Fluxes(domain::Domain)

  # Get parameters.
  (; nx, ny, nz) = domain

  # Initialize the fluxes.
  (phirho, phirhop, phiu, phiv, phiw) = (
    OffsetArray(
      zeros((nx + 2, ny + 2, nz + 2, 3)),
      (-1):(nx),
      (-1):(ny),
      (-1):(nz),
      1:3,
    ) for i in 1:5
  )

  # Return a Fluxes instance.
  return Fluxes(phirho, phirhop, phiu, phiv, phiw)
end
