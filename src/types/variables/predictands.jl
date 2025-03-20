struct Predictands{A <: AbstractArray{<:AbstractFloat, 3}}
  rho::A
  rhop::A
  u::A
  v::A
  w::A
  pip::A
end

function Predictands(
  namelists::Namelists,
  constants::Constants,
  domain::Domain,
  model::PseudoIncompressible,
  testcase::MountainWave,
)

  # Get parameters.
  (; backgroundflow_dim) = namelists.atmosphere
  (; nbx, nby, nbz) = namelists.domain
  (; uref) = constants
  (; nx, ny, nz) = domain

  # Initialize the predictands.
  (rho, rhop, u, v, w, pip) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:6
  )

  # Initial winds.
  u .= backgroundflow_dim[1] / uref
  v .= backgroundflow_dim[2] / uref
  w .= backgroundflow_dim[3] / uref

  # Return a Predictands instance.
  return Predictands(rho, rhop, u, v, w, pip)
end
