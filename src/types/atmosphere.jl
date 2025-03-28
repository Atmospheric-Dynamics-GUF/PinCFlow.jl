struct Atmosphere{
  A <: AbstractArray{<:AbstractFloat, 3},
  B <: AbstractVector{<:AbstractFloat},
}
  pstrattfc::A
  thetastrattfc::A
  rhostrattfc::A
  bvsstrattfc::A
  f_cor_nd::B
end

function Atmosphere(
  namelists::Namelists,
  constants::Constants,
  domain::Domain,
  grid::Grid,
)
  (; model) = namelists.setting
  (; background, corset) = namelists.atmosphere
  return Atmosphere(
    namelists,
    constants,
    domain,
    grid,
    model,
    background,
    corset,
  )
end

function Atmosphere(
  namelists::Namelists,
  constants::Constants,
  domain::Domain,
  grid::Grid,
  model::PseudoIncompressible,
  background::Isothermal,
  corset::ConstantCoriolis,
)
  # Get parameters.
  (; nbx, nby) = namelists.domain
  (; temp0_dim, press0_dim, f_coriolis_dim) = namelists.atmosphere
  (; tref, thetaref, pref, kappa, sig, gamma, g_ndim) = constants
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

  # Set Coriolis parameter.
  f_cor_nd = OffsetArray(zeros(ny + 2), 0:(ny + 1))
  f_cor_nd[0:(ny + 1)] .= f_coriolis_dim .* tref

  # Return an Atmosphere instance.
  return Atmosphere(
    pstrattfc,
    thetastrattfc,
    rhostrattfc,
    bvsstrattfc,
    f_cor_nd,
  )
end
