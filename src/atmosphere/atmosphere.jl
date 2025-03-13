using SimpleUnPack
using OffsetArrays

struct Atmosphere{A <: OffsetArray{<:AbstractFloat, 3}, B <: AbstractFloat}

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


function map(level, lz)
    # Vertical grid stretching.
    stretch_exponent = 1.0 # TODO
    if level < 0
        return -lz[1] * (-level / lz[1])^stretch_exponent
    elseif level > lz[1]
        return 2 * lz[1] - lz[1] * ((2 * lz[1] - level) / lz[1])^stretch_exponent
    else
        return lz[1] * (level / lz[1])^stretch_exponent
    end
end

function vertWind(i, j, k, model)
    # Transformation of the vertical wind.
    (; u, v, w) = model.variables.prognostic_fields

    uEdgeR = u[i, j, k]
    uUEdgeR = u[i, j, k + 1]
    uEdgeL = u[i - 1, j, k]
    uUEdgeL = u[i - 1, j, k + 1]
    vEdgeF = v[i, j, k]
    vUEdgeF = v[i, j, k + 1]
    vEdgeB = v[i, j - 1, k]
    vUEdgeB = v[i, j - 1, k + 1]
    wEdgeU = w[i, j, k]

    return trafo(i, j, k, uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, vUEdgeF, vEdgeB,
                 vUEdgeB,
                 wEdgeU, "car", model)
end

function trafo(i, j, k, uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, vUEdgeF, vEdgeB, vUEdgeB,
               wEdgeU, wind, model)
    # Assuming jac and met are defined elsewhere
    # Define variables as in the original code

    (; jac, met) = model.grid

    jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) + jac(i, j, k + 1))
    uC = 0.5 * (uEdgeR + uEdgeL)
    uU = 0.5 * (uUEdgeR + uUEdgeL)
    vC = 0.5 * (vEdgeF + vEdgeB)
    vU = 0.5 * (vUEdgeF + vUEdgeB)

    if wind == "car"
        trafoTFC = jacEdgeU * (-(jac(i, j, k + 1) *
                      (met(i, j, k, 1, 3) * uC + met(i, j, k, 2, 3) * vC) +
                      jac(i, j, k) *
                      (met(i, j, k + 1, 1, 3) * uU + met(i, j, k + 1, 2, 3) * vU)) /
                    (jac(i, j, k) + jac(i, j, k + 1)) + wEdgeU)
    elseif wind == "tfc"
        trafoTFC = (jac(i, j, k + 1) * (met(i, j, k, 1, 3) * uC + met(i, j, k, 2, 3) * vC) +
                    jac(i, j, k) *
                    (met(i, j, k + 1, 1, 3) * uU + met(i, j, k + 1, 2, 3) * vU)) /
                   (jac(i, j, k) + jac(i, j, k + 1)) + wEdgeU / jacEdgeU
    end
    return trafoTFC
end

function stressTensTFC(i, j, k, mu, nu, semi)
    (; var, cache) = semi
    (; jac) = cache
    # Assuming jac, met, and vertWind are functions defined elsewhere

    # Define variables as in the original code
    jacEdgeR = 0.5 * (jac(i, j, k) + jac(i + 1, j, k))
    jacEdgeL = 0.5 * (jac(i, j, k) + jac(i - 1, j, k))
    jacEdgeF = 0.5 * (jac(i, j, k) + jac(i, j + 1, k))
    jacEdgeB = 0.5 * (jac(i, j, k) + jac(i, j - 1, k))
    jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) + jac(i, j, k + 1))
    jacEdgeD = 2.0 * jac(i, j, k) * jac(i, j, k - 1) / (jac(i, j, k) + jac(i, j, k - 1))

    # Accessing array elements in var
    uF = 0.5 * (var.u[i, j + 1, k] + var.u[i - 1, j + 1, k])
    uB = 0.5 * (var.u[i, j - 1, k] + var.u[i - 1, j - 1, k])
    uU = 0.5 * (var.u[i, j, k + 1] + var.u[i - 1, j, k + 1])
    uD = 0.5 * (var.u[i, j, k - 1] + var.u[i - 1, j, k - 1])

    vR = 0.5 * (var.v[i + 1, j, k] + var.v[i + 1, j - 1, k])
    vL = 0.5 * (var.v[i - 1, j, k] + var.v[i - 1, j - 1, k])
    vU = 0.5 * (var.v[i, j, k + 1] + var.v[i, j - 1, k + 1])
    vD = 0.5 * (var.v[i, j, k - 1] + var.v[i, j - 1, k - 1])

    wR = 0.5 * (vertWind(i + 1, j, k, semi) + vertWind(i + 1, j, k - 1, semi))
    wL = 0.5 * (vertWind(i - 1, j, k, semi) + vertWind(i - 1, j, k - 1, semi))
    wF = 0.5 * (vertWind(i, j + 1, k, semi) + vertWind(i, j + 1, k - 1, semi))
    wB = 0.5 * (vertWind(i, j - 1, k, semi) + vertWind(i, j - 1, k - 1, semi))

    # Conditional logic for stress tensor calculation
    if mu == 1 && nu == 1
        stressTensTFC = 2.0 * (var.u[i, j, k] - var.u[i - 1, j, k]) / dx +
                        met(i, j, k, 1, 3) * (uU - uD) / dz -
                        2.0 / 3.0 *
                        ((jacEdgeR * var.u[i, j, k] - jacEdgeL * var.u[i - 1, j, k]) / dx +
                         (jacEdgeF * var.v[i, j, k] - jacEdgeB * var.v[i, j - 1, k]) / dy +
                         (jacEdgeU * var.w[i, j, k] - jacEdgeD * var.w[i, j, k - 1]) / dz) /
                        jac(i, j, k)
    elseif (mu == 1 && nu == 2) || (mu == 2 && nu == 1)
        stressTensTFC = 0.5 * (uF - uB) / dy +
                        0.5 * met(i, j, k, 2, 3) * (uU - uD) / dz +
                        0.5 * (vR - vL) / dx +
                        0.5 * met(i, j, k, 1, 3) * (vU - vD) / dz
    elseif (mu == 1 && nu == 3) || (mu == 3 && nu == 1)
        stressTensTFC = 0.5 * (uU - uD) / dz / jac(i, j, k) +
                        0.5 * (wR - wL) / dx +
                        met(i, j, k, 1, 3) *
                        (vertWind(i, j, k, semi) - vertWind(i, j, k - 1, semi)) /
                        dz
    elseif mu == 2 && nu == 2
        stressTensTFC = 2.0 * (var.v[i, j, k] - var.v[i, j - 1, k]) / dy +
                        met(i, j, k, 2, 3) * (vU - vD) / dz -
                        2.0 / 3.0 *
                        ((jacEdgeR * var.u[i, j, k] - jacEdgeL * var.u[i - 1, j, k]) / dx +
                         (jacEdgeF * var.v[i, j, k] - jacEdgeB * var.v[i, j - 1, k]) / dy +
                         (jacEdgeU * var.w[i, j, k] - jacEdgeD * var.w[i, j, k - 1]) / dz) /
                        jac(i, j, k)
    elseif (mu == 2 && nu == 3) || (mu == 3 && nu == 2)
        stressTensTFC = (0.5 * (vU - vD) / dz / jac(i, j, k) +
                         0.5 * (wF - wB) / dy +
                         met(i, j, k, 2, 3) *
                         (vertWind(i, j, k, semi) - vertWind(i, j, k - 1, semi)) / dz)
    elseif mu == 3 && nu == 3
        stressTensTFC = 2.0 * (vertWind(i, j, k, semi) - vertWind(i, j, k - 1, semi)) / dz /
                        jac(i, j, k) -
                        2.0 / 3.0 *
                        ((jacEdgeR * var.u[i, j, k] - jacEdgeL * var.u[i - 1, j, k]) / dx +
                         (jacEdgeF * var.v[i, j, k] - jacEdgeB * var.v[i, j - 1, k]) / dy +
                         (jacEdgeU * var.w[i, j, k] - jacEdgeD * var.w[i, j, k - 1]) / dz) /
                        jac(i, j, k)
    end

    return stressTensTFC
end

function set_topography_boundary!(model, field)
    (; nx, ny) = model.domain
    (; nbx, nby) = model.parameters.domain

    for i in 1:nbx
        field[nx + i, :] = field[i, :]
        field[-i + 1, :] = field[nx - i + 1, :]
    end
    for j in 1:nby
        field[:, ny + j] = field[:, j]
        field[:, -j + 1] = field[:, ny - j + 1]
    end
end
