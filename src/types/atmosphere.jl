struct Atmosphere{A <: AbstractArray{<:AbstractFloat, 3}}
  pstrattfc::A
  thetastrattfc::A
  rhostrattfc::A
  bvsstrattfc::A
end

function Atmosphere(
  namelists::Namelists,
  constants::Constants,
  domain::Domain,
  grid::Grid,
  model::PseudoIncompressible,
  background::Isothermal,
)
  # Get parameters.
  (; nbx, nby) = namelists.domain
  (; temp0_dim, press0_dim) = namelists.atmosphere
  (; thetaref, pref, ma, fr, kappa, sig, gamma, g_ndim) = constants
  (; nx, ny, nz) = domain
  (; ztfc, jac, dz) = grid

  # Initialize the background fields.
  (pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 4)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-1):(nz + 2),
    ) for i in 1:4
  )

  t0 = temp0_dim / thetaref
  p0 = press0_dim / pref

  # Define 3D background fields.
  for k in (-1):(nz + 2)
    # Define pStratTFC.
    pstrattfc[:, :, k] .= p0 .* exp.(-sig .* ztfc[:, :, k] ./ gamma ./ t0)
    # Define thetaStratTFC.
    thetastrattfc[:, :, k] .= t0 .* exp.(kappa .* sig ./ t0 .* ztfc[:, :, k])
    # Define rhoStratTFC.
    rhostrattfc[:, :, k] .= pstrattfc[:, :, k] ./ thetastrattfc[:, :, k]
  end

  # Define bvsStratTFC.
  bvsstrattfc .= 0.0
  # Lower boundary.
  bvsstrattfc[:, :, -1] .=
    g_ndim ./ thetastrattfc[:, :, 0] ./ jac[:, :, 0] .*
    (thetastrattfc[:, :, 1] .- thetastrattfc[:, :, 0]) ./ dz
  bvsstrattfc[:, :, 0] .= bvsstrattfc[:, :, -1]
  # Between boundaries.
  for k in 1:nz
    bvsstrattfc[:, :, k] .=
      g_ndim ./ thetastrattfc[:, :, k] ./ jac[:, :, k] .* 0.5 .*
      (thetastrattfc[:, :, k + 1] .- thetastrattfc[:, :, k - 1]) ./ dz
  end
  # Upper boundary.
  bvsstrattfc[:, :, nz + 1] .=
    g_ndim ./ thetastrattfc[:, :, nz + 1] ./ jac[:, :, nz + 1] .*
    (thetastrattfc[:, :, nz + 1] .- thetastrattfc[:, :, nz]) ./ dz
  bvsstrattfc[:, :, nz + 2] .= bvsstrattfc[:, :, nz + 1]

  # Return an Atmosphere instance.
  return Atmosphere(
    pstrattfc,
    thetastrattfc,
    rhostrattfc,
    bvsstrattfc,
  )
end
