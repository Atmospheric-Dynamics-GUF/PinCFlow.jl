struct Auxiliaries{
  A <: AbstractArray{<:AbstractFloat, 3},
  B <: AbstractVector{<:AbstractFloat},
  C <: AbstractVector{<:AbstractFloat},
  D <: AbstractVector{<:AbstractFloat},
  E <: AbstractMatrix{<:AbstractFloat},
  F <: AbstractMatrix{<:AbstractFloat},
  G <: AbstractMatrix{<:AbstractFloat},
}
  phi::A
  phix::B
  phiy::C
  phiz::D
  phitildex::E
  phitildey::F
  phitildez::G
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
  phix = zeros(nx + 2 * nbx + 1)
  phiy = zeros(ny + 2 * nby + 1)
  phiz = zeros(nz + 2 * nbz + 1)
  phitildex = zeros(nx + 2 * nbx + 1, 2)
  phitildey = zeros(ny + 2 * nby + 1, 2)
  phitildez = zeros(nz + 2 * nbz + 1, 2)

  # Return an Auxiliaries instance.
  return Auxiliaries(
    phi,
    phix,
    phiy,
    phiz,
    phitildex,
    phitildey,
    phitildez,
  )
end
