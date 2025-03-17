struct Reconstructions{A <: OffsetArray{<:AbstractFloat, 5}}
  rhotilde::A
  rhoptilde::A
  utilde::A
  vtilde::A
  wtilde::A
end

function Reconstructions(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the reconstructed variables.
  (rhotilde, rhoptilde, utilde, vtilde, wtilde) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1, 3, 2)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
      1:3,
      0:1,
    ) for i in 1:5
  )

  # Return a Reconstructions instance.
  return Reconstructions(rhotilde, rhoptilde, utilde, vtilde, wtilde)
end
