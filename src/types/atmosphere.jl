struct Atmosphere{A <: AbstractArray{<:AbstractFloat, 3}, B <: AbstractFloat}

  # Reference atmosphere.
  pstrattfc::A
  thetastrattfc::A
  rhostrattfc::A
  bvsstrattfc::A

  # Scaled reference values.
  n2::B
  nn::B
  t0::B
  p0::B
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
  n2 = ma^2 / fr^4 * kappa / t0
  nn = sqrt(n2)

  # Define 3D background fields.
  for i in (-nbx):(nx + nbx)
    for j in (-nby):(ny + nby)
      for k in (-1):(nz + 2)
        # Define pStratTFC.
        pstrattfc[i, j, k] = p0 * exp(-sig * ztfc[i, j, k] / gamma / t0)
        # Define thetaStratTFC.
        thetastrattfc[i, j, k] = t0 * exp(kappa * sig / t0 * ztfc[i, j, k])
        # Define rhoStratTFC.
        rhostrattfc[i, j, k] = pstrattfc[i, j, k] / thetastrattfc[i, j, k]
      end
    end
  end

  # Define bvsStratTFC.
  bvsstrattfc .= 0.0
  for i in (-nbx):(nx + nbx)
    for j in (-nby):(ny + nby)
      # Lower boundary.
      bvsstrattfc[i, j, -1] =
        g_ndim / thetastrattfc[i, j, 0] / jac[i, j, 0] *
        (thetastrattfc[i, j, 1] - thetastrattfc[i, j, 0]) / dz
      bvsstrattfc[i, j, 0] = bvsstrattfc[i, j, -1]
      # Between boundaries.
      for k in 1:(nz)
        bvsstrattfc[i, j, k] =
          g_ndim / thetastrattfc[i, j, k] / jac[i, j, k] *
          0.5 *
          (thetastrattfc[i, j, k + 1] - thetastrattfc[i, j, k - 1]) / dz
      end
      # Upper boundary.
      bvsstrattfc[i, j, nz + 1] =
        g_ndim / thetastrattfc[i, j, nz + 1] / jac[i, j, nz + 1] *
        (thetastrattfc[i, j, nz + 1] - thetastrattfc[i, j, nz]) / dz
      bvsstrattfc[i, j, nz + 2] = bvsstrattfc[i, j, nz + 1]
    end
  end

  # Return an Atmosphere instance.
  return Atmosphere(
    pstrattfc,
    thetastrattfc,
    rhostrattfc,
    bvsstrattfc,
    n2,
    nn,
    t0,
    p0,
  )
end
