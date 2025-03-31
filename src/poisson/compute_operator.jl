function compute_operator!(
  state::State,
  dt::AbstractFloat,
  opt::EXPL,
  facray::AbstractFloat,
)
  (; preconditioner) = state.namelists.poisson
  (; zboundaries) = state.namelists.setting
  (; i0, i1, j0, j1, k0, k1) = state.domain
  (; dx, dy, dz, jac, met) = state.grid
  (; pstrattfc, rhostrattfc) = state.atmosphere
  (;
    ac_b,
    acv_b,
    ach_b,
    al_b,
    ar_b,
    ab_b,
    af_b,
    ad_b,
    au_b,
    aru_b,
    ard_b,
    alu_b,
    ald_b,
    afu_b,
    afd_b,
    abu_b,
    abd_b,
    auu_b,
    add_b,
    aruu_b,
    ardd_b,
    aluu_b,
    aldd_b,
    afuu_b,
    afdd_b,
    abuu_b,
    abdd_b,
  ) = state.poisson.tensor
  (; rho) = state.variables.predictands

  # Compute tensor elements for TFC.
  for k in k0:k1, j in j0:j1, i in i0:i1
    # Compute scaling factors.
    fcscal = sqrt(pstrattfc[i, j, k]^2.0 / rhostrattfc[i, j, k])
    fcscal_r = sqrt(pstrattfc[i + 1, j, k]^2.0 / rhostrattfc[i + 1, j, k])
    fcscal_l = sqrt(pstrattfc[i - 1, j, k]^2.0 / rhostrattfc[i - 1, j, k])
    fcscal_f = sqrt(pstrattfc[i, j + 1, k]^2.0 / rhostrattfc[i, j + 1, k])
    fcscal_b = sqrt(pstrattfc[i, j - 1, k]^2.0 / rhostrattfc[i, j - 1, k])
    fcscal_u = sqrt(pstrattfc[i, j, k + 1]^2.0 / rhostrattfc[i, j, k + 1])
    fcscal_d = sqrt(pstrattfc[i, j, k - 1]^2.0 / rhostrattfc[i, j, k - 1])
    fcscal_ru =
      sqrt(pstrattfc[i + 1, j, k + 1]^2.0 / rhostrattfc[i + 1, j, k + 1])
    fcscal_rd =
      sqrt(pstrattfc[i + 1, j, k - 1]^2.0 / rhostrattfc[i + 1, j, k - 1])
    fcscal_lu =
      sqrt(pstrattfc[i - 1, j, k + 1]^2.0 / rhostrattfc[i - 1, j, k + 1])
    fcscal_ld =
      sqrt(pstrattfc[i - 1, j, k - 1]^2.0 / rhostrattfc[i - 1, j, k - 1])
    fcscal_fu =
      sqrt(pstrattfc[i, j + 1, k + 1]^2.0 / rhostrattfc[i, j + 1, k + 1])
    fcscal_fd =
      sqrt(pstrattfc[i, j + 1, k - 1]^2.0 / rhostrattfc[i, j + 1, k - 1])
    fcscal_bu =
      sqrt(pstrattfc[i, j - 1, k + 1]^2.0 / rhostrattfc[i, j - 1, k + 1])
    fcscal_bd =
      sqrt(pstrattfc[i, j - 1, k - 1]^2.0 / rhostrattfc[i, j - 1, k - 1])
    fcscal_uu = sqrt(pstrattfc[i, j, k + 2]^2.0 / rhostrattfc[i, j, k + 2])
    fcscal_dd = sqrt(pstrattfc[i, j, k - 2]^2.0 / rhostrattfc[i, j, k - 2])
    fcscal_ruu =
      sqrt(pstrattfc[i + 1, j, k + 2]^2.0 / rhostrattfc[i + 1, j, k + 2])
    fcscal_rdd =
      sqrt(pstrattfc[i + 1, j, k - 2]^2.0 / rhostrattfc[i + 1, j, k - 2])
    fcscal_luu =
      sqrt(pstrattfc[i - 1, j, k + 2]^2.0 / rhostrattfc[i - 1, j, k + 2])
    fcscal_ldd =
      sqrt(pstrattfc[i - 1, j, k - 2]^2.0 / rhostrattfc[i - 1, j, k - 2])
    fcscal_fuu =
      sqrt(pstrattfc[i, j + 1, k + 2]^2.0 / rhostrattfc[i, j + 1, k + 2])
    fcscal_fdd =
      sqrt(pstrattfc[i, j + 1, k - 2]^2.0 / rhostrattfc[i, j + 1, k - 2])
    fcscal_buu =
      sqrt(pstrattfc[i, j - 1, k + 2]^2.0 / rhostrattfc[i, j - 1, k + 2])
    fcscal_bdd =
      sqrt(pstrattfc[i, j - 1, k - 2]^2.0 / rhostrattfc[i, j - 1, k - 2])

    # Compute inverse Jacobian.
    jacinv = 1.0 / jac[i, j, k]

    # Compute P coefficients (divergence).
    pedgerdiv =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
      )
    pedgeldiv =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i - 1, j, k] * pstrattfc[i - 1, j, k]
      )
    pedgefdiv =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
      )
    pedgebdiv =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j - 1, k] * pstrattfc[i, j - 1, k]
      )
    pedgeudiv =
      jac[i, j, k] *
      jac[i, j, k + 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1])
    pedgeddiv =
      jac[i, j, k] *
      jac[i, j, k - 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k - 1]) /
      (jac[i, j, k] + jac[i, j, k - 1])

    # Compute P coefficients (pressure gradient).
    pedgergra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
    pedgelgra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i - 1, j, k])
    pedgefgra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
    pedgebgra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j - 1, k])
    pedgeugra =
      (
        jac[i, j, k + 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    pedgedgra =
      (
        jac[i, j, k - 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k - 1]
      ) / (jac[i, j, k] + jac[i, j, k - 1])

    # Compute density coefficients.
    rhoedger =
      0.5 * (
        rho[i, j, k] +
        rho[i + 1, j, k] +
        rhostrattfc[i, j, k] +
        rhostrattfc[i + 1, j, k]
      )
    rhoedgel =
      0.5 * (
        rho[i, j, k] +
        rho[i - 1, j, k] +
        rhostrattfc[i, j, k] +
        rhostrattfc[i - 1, j, k]
      )
    rhoedgef =
      0.5 * (
        rho[i, j, k] +
        rho[i, j + 1, k] +
        rhostrattfc[i, j, k] +
        rhostrattfc[i, j + 1, k]
      )
    rhoedgeb =
      0.5 * (
        rho[i, j, k] +
        rho[i, j - 1, k] +
        rhostrattfc[i, j, k] +
        rhostrattfc[i, j - 1, k]
      )
    rhoedgeu =
      (
        jac[i, j, k + 1] * (rho[i, j, k] + rhostrattfc[i, j, k]) +
        jac[i, j, k] * (rho[i, j, k + 1] + rhostrattfc[i, j, k + 1])
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    rhoedged =
      (
        jac[i, j, k - 1] * (rho[i, j, k] + rhostrattfc[i, j, k]) +
        jac[i, j, k] * (rho[i, j, k - 1] + rhostrattfc[i, j, k - 1])
      ) / (jac[i, j, k] + jac[i, j, k - 1])

    # Interpolate metric-tensor elements.
    met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])
    met13edgel = 0.5 * (met[i, j, k, 1, 3] + met[i - 1, j, k, 1, 3])
    met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])
    met23edgeb = 0.5 * (met[i, j, k, 2, 3] + met[i, j - 1, k, 2, 3])
    met13edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 1, 3] +
        jac[i, j, k] * met[i, j, k + 1, 1, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met23edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 2, 3] +
        jac[i, j, k] * met[i, j, k + 1, 2, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met33edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 3, 3] +
        jac[i, j, k] * met[i, j, k + 1, 3, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met13edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 1, 3] +
        jac[i, j, k] * met[i, j, k - 1, 1, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    met23edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 2, 3] +
        jac[i, j, k] * met[i, j, k - 1, 2, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    met33edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 3, 3] +
        jac[i, j, k] * met[i, j, k - 1, 3, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])

    # --------------------- A(i,j,k) ---------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      ac =
        -jacinv / dx * (
          pedgerdiv / rhoedger *
          pedgergra *
          (1.0 / dx + 0.75 * met13edger / dz) +
          pedgeldiv / rhoedgel *
          pedgelgra *
          (1.0 / dx - 0.75 * met13edgel / dz)
        ) -
        jacinv / dy * (
          pedgefdiv / rhoedgef *
          pedgefgra *
          (1.0 / dy + 0.75 * met23edgef / dz) +
          pedgebdiv / rhoedgeb *
          pedgebgra *
          (1.0 / dy - 0.75 * met23edgeb / dz)
        ) - jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met33edgeu / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ac =
        -jacinv / dx * (
          pedgerdiv / rhoedger *
          pedgergra *
          (1.0 / dx - 0.75 * met13edger / dz) +
          pedgeldiv / rhoedgel *
          pedgelgra *
          (1.0 / dx + 0.75 * met13edgel / dz)
        ) -
        jacinv / dy * (
          pedgefdiv / rhoedgef *
          pedgefgra *
          (1.0 / dy - 0.75 * met23edgef / dz) +
          pedgebdiv / rhoedgeb *
          pedgebgra *
          (1.0 / dy + 0.75 * met23edgeb / dz)
        ) - jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met33edged / dz
    else
      ac =
        -jacinv / dx * (
          pedgerdiv / rhoedger * pedgergra / dx +
          pedgeldiv / rhoedgel * pedgelgra / dx
        ) -
        jacinv / dy * (
          pedgefdiv / rhoedgef * pedgefgra / dy +
          pedgebdiv / rhoedgeb * pedgebgra / dy
        ) -
        jacinv / dz * (
          pedgeudiv / rhoedgeu * pedgeugra * met33edgeu / dz +
          pedgeddiv / rhoedged * pedgedgra * met33edged / dz
        )
    end
    # -------------------- A(i+1,j,k) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      ar =
        jacinv / dx * pedgerdiv / rhoedger *
        pedgergra *
        (1.0 / dx - 0.75 * met13edger / dz) +
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * 0.5 * met13edgeu / dx *
        jac[i + 1, j, k + 1] / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ar =
        jacinv / dx * pedgerdiv / rhoedger *
        pedgergra *
        (1.0 / dx + 0.75 * met13edger / dz) -
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * 0.5 * met13edged / dx *
        jac[i + 1, j, k - 1] / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
    else
      ar =
        jacinv / dx * pedgerdiv / rhoedger * pedgergra / dx +
        jacinv / dz * (
          pedgeudiv / rhoedgeu * pedgeugra * met13edgeu * 0.5 / dx *
          jac[i + 1, j, k + 1] / (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) -
          pedgeddiv / rhoedged * pedgedgra * met13edged * 0.5 / dx *
          jac[i + 1, j, k - 1] / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
        )
    end

    # -------------------- A(i-1,j,k) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      al =
        jacinv / dx * pedgeldiv / rhoedgel *
        pedgelgra *
        (1.0 / dx + 0.75 * met13edgel / dz) -
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * 0.5 * met13edgeu / dx *
        jac[i - 1, j, k + 1] / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      al =
        jacinv / dx * pedgeldiv / rhoedgel *
        pedgelgra *
        (1.0 / dx - 0.75 * met13edgel / dz) +
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * 0.5 * met13edged / dx *
        jac[i - 1, j, k - 1] / (jac[i - 1, j, k] + jac[i - 1, j, k - 1])
    else
      al =
        jacinv / dx * pedgeldiv / rhoedgel * pedgelgra / dx -
        jacinv / dz * (
          pedgeudiv / rhoedgeu * pedgeugra * met13edgeu * 0.5 / dx *
          jac[i - 1, j, k + 1] / (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
          pedgeddiv / rhoedged * pedgedgra * met13edged * 0.5 / dx *
          jac[i - 1, j, k - 1] / (jac[i - 1, j, k] + jac[i - 1, j, k - 1])
        )
    end

    # -------------------- A(i,j+1,k) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      af =
        jacinv / dy * pedgefdiv / rhoedgef *
        pedgefgra *
        (1.0 / dy - 0.75 * met23edgef / dz) +
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * 0.5 * met23edgeu / dy *
        jac[i, j + 1, k + 1] / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      af =
        jacinv / dy * pedgefdiv / rhoedgef *
        pedgefgra *
        (1.0 / dy + 0.75 * met23edgef / dz) -
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * 0.5 * met23edged / dy *
        jac[i, j + 1, k - 1] / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
    else
      af =
        jacinv / dy * pedgefdiv / rhoedgef * pedgefgra / dy +
        jacinv / dz * (
          pedgeudiv / rhoedgeu * pedgeugra * met23edgeu * 0.5 / dy *
          jac[i, j + 1, k + 1] / (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) -
          pedgeddiv / rhoedged * pedgedgra * met23edged * 0.5 / dy *
          jac[i, j + 1, k - 1] / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
        )
    end

    # -------------------- A(i,j-1,k) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      ab =
        jacinv / dy * pedgebdiv / rhoedgeb *
        pedgebgra *
        (1.0 / dy + 0.75 * met23edgeb / dz) -
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * 0.5 * met23edgeu / dy *
        jac[i, j - 1, k + 1] / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ab =
        jacinv / dy * pedgebdiv / rhoedgeb *
        pedgebgra *
        (1.0 / dy - 0.75 * met23edgeb / dz) +
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * 0.5 * met23edged / dy *
        jac[i, j - 1, k - 1] / (jac[i, j - 1, k] + jac[i, j - 1, k - 1])
    else
      ab =
        jacinv / dy * pedgebdiv / rhoedgeb * pedgebgra / dy -
        jacinv / dz * (
          pedgeudiv / rhoedgeu * pedgeugra * met23edgeu * 0.5 / dy *
          jac[i, j - 1, k + 1] / (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
          pedgeddiv / rhoedged * pedgedgra * met23edged * 0.5 / dy *
          jac[i, j - 1, k - 1] / (jac[i, j - 1, k] + jac[i, j - 1, k - 1])
        )
    end

    # -------------------- A(i,j,k+1) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      au =
        jacinv / dx * (
          pedgerdiv / rhoedger * pedgergra * met13edger / dz -
          pedgeldiv / rhoedgel * pedgelgra * met13edgel / dz
        ) +
        jacinv / dy * (
          pedgefdiv / rhoedgef * pedgefgra * met23edgef / dz -
          pedgebdiv / rhoedgeb * pedgebgra * met23edgeb / dz
        ) +
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met33edgeu / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      au = 0.0
    else
      au =
        jacinv / dx * (
          pedgerdiv / rhoedger * pedgergra * met13edger * 0.25 / dz -
          pedgeldiv / rhoedgel * pedgelgra * met13edgel * 0.25 / dz
        ) +
        jacinv / dy * (
          pedgefdiv / rhoedgef * pedgefgra * met23edgef * 0.25 / dz -
          pedgebdiv / rhoedgeb * pedgebgra * met23edgeb * 0.25 / dz
        ) +
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met33edgeu / dz
    end

    # -------------------- A(i,j,k-1) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      ad = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ad =
        -jacinv / dx * (
          pedgerdiv / rhoedger * pedgergra * met13edger / dz -
          pedgeldiv / rhoedgel * pedgelgra * met13edgel / dz
        ) -
        jacinv / dy * (
          pedgefdiv / rhoedgef * pedgefgra * met23edgef / dz -
          pedgebdiv / rhoedgeb * pedgebgra * met23edgeb / dz
        ) + jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met33edged / dz
    else
      ad =
        -jacinv / dx * (
          pedgerdiv / rhoedger * pedgergra * met13edger * 0.25 / dz -
          pedgeldiv / rhoedgel * pedgelgra * met13edgel * 0.25 / dz
        ) -
        jacinv / dy * (
          pedgefdiv / rhoedgef * pedgefgra * met23edgef * 0.25 / dz -
          pedgebdiv / rhoedgeb * pedgebgra * met23edgeb * 0.25 / dz
        ) + jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met33edged / dz
    end

    # ------------------- A(i+1,j,k+1) -------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      aru =
        jacinv / dx * pedgerdiv / rhoedger * pedgergra * met13edger / dz +
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met13edgeu * 0.5 / dx *
        jac[i + 1, j, k] / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      aru = 0.0
    else
      aru =
        jacinv / dx * pedgerdiv / rhoedger * pedgergra * met13edger * 0.25 /
        dz +
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met13edgeu * 0.5 / dx *
        jac[i + 1, j, k] / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
    end

    # ------------------- A(i+1,j,k-1) -------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      ard = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ard =
        -jacinv / dx * pedgerdiv / rhoedger * pedgergra * met13edger / dz -
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met13edged * 0.5 / dx *
        jac[i + 1, j, k] / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
    else
      ard =
        -jacinv / dx * pedgerdiv / rhoedger * pedgergra * met13edger * 0.25 /
        dz -
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met13edged * 0.5 / dx *
        jac[i + 1, j, k] / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
    end

    # ------------------- A(i-1,j,k+1) -------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      alu =
        -jacinv / dx * pedgeldiv / rhoedgel * pedgelgra * met13edgel / dz -
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met13edgeu * 0.5 / dx *
        jac[i - 1, j, k] / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      alu = 0.0
    else
      alu =
        -jacinv / dx * pedgeldiv / rhoedgel * pedgelgra * met13edgel * 0.25 /
        dz -
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met13edgeu * 0.5 / dx *
        jac[i - 1, j, k] / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
    end

    # ------------------- A(i-1,j,k-1) -------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      ald = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ald =
        jacinv / dx * pedgeldiv / rhoedgel * pedgelgra * met13edgel / dz +
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met13edged * 0.5 / dx *
        jac[i - 1, j, k] / (jac[i - 1, j, k] + jac[i - 1, j, k - 1])
    else
      ald =
        jacinv / dx * pedgeldiv / rhoedgel * pedgelgra * met13edgel * 0.25 /
        dz +
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met13edged * 0.5 / dx *
        jac[i - 1, j, k] / (jac[i - 1, j, k] + jac[i - 1, j, k - 1])
    end

    # ------------------- A(i,j+1,k+1) -------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      afu =
        jacinv / dy * pedgefdiv / rhoedgef * pedgefgra * met23edgef / dz +
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met23edgeu * 0.5 / dy *
        jac[i, j + 1, k] / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      afu = 0.0
    else
      afu =
        jacinv / dy * pedgefdiv / rhoedgef * pedgefgra * met23edgef * 0.25 /
        dz +
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met23edgeu * 0.5 / dy *
        jac[i, j + 1, k] / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
    end

    # ------------------- A(i,j+1,k-1) -------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      afd = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      afd =
        -jacinv / dy * pedgefdiv / rhoedgef * pedgefgra * met23edgef / dz -
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met23edged * 0.5 / dy *
        jac[i, j + 1, k] / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
    else
      afd =
        -jacinv / dy * pedgefdiv / rhoedgef * pedgefgra * met23edgef * 0.25 /
        dz -
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met23edged * 0.5 / dy *
        jac[i, j + 1, k] / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
    end

    # ------------------- A(i,j-1,k+1) -------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      abu =
        -jacinv / dy * pedgebdiv / rhoedgeb * pedgebgra * met23edgeb / dz -
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met23edgeu * 0.5 / dy *
        jac[i, j - 1, k] / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      abu = 0.0
    else
      abu =
        -jacinv / dy * pedgebdiv / rhoedgeb * pedgebgra * met23edgeb * 0.25 /
        dz -
        jacinv / dz * pedgeudiv / rhoedgeu * pedgeugra * met23edgeu * 0.5 / dy *
        jac[i, j - 1, k] / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])
    end

    # ------------------- A(i,j-1,k-1) -------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      abd = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      abd =
        jacinv / dy * pedgebdiv / rhoedgeb * pedgebgra * met23edgeb / dz +
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met23edged * 0.5 / dy *
        jac[i, j - 1, k] / (jac[i, j - 1, k] + jac[i, j - 1, k - 1])
    else
      abd =
        jacinv / dy * pedgebdiv / rhoedgeb * pedgebgra * met23edgeb * 0.25 /
        dz +
        jacinv / dz * pedgeddiv / rhoedged * pedgedgra * met23edged * 0.5 / dy *
        jac[i, j - 1, k] / (jac[i, j - 1, k] + jac[i, j - 1, k - 1])
    end

    # ------------------- A(i,j,k+2) ---------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      auu =
        -jacinv / dx * (
          pedgerdiv / rhoedger * pedgergra * 0.25 * met13edger / dz -
          pedgeldiv / rhoedgel * pedgelgra * 0.25 * met13edgel / dz
        ) -
        jacinv / dy * (
          pedgefdiv / rhoedgef * pedgefgra * 0.25 * met23edgef / dz -
          pedgebdiv / rhoedgeb * pedgebgra * 0.25 * met23edgeb / dz
        )
    else
      auu = 0.0
    end

    # ------------------- A(i,j,k-2) ---------------------

    if k == k1 && zboundaries == SolidWallBoundaries()
      add =
        jacinv / dx * (
          pedgerdiv / rhoedger * pedgergra * 0.25 * met13edger / dz -
          pedgeldiv / rhoedgel * pedgelgra * 0.25 * met13edgel / dz
        ) +
        jacinv / dy * (
          pedgefdiv / rhoedgef * pedgefgra * 0.25 * met23edgef / dz -
          pedgebdiv / rhoedgeb * pedgebgra * 0.25 * met23edgeb / dz
        )
    else
      add = 0.0
    end

    # ------------------ A(i+1,j,k+2) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      aruu =
        -jacinv / dx * pedgerdiv / rhoedger * pedgergra * 0.25 * met13edger / dz
    else
      aruu = 0.0
    end

    # ------------------ A(i+1,j,k-2) --------------------

    if k == k1 && zboundaries == SolidWallBoundaries()
      ardd =
        jacinv / dx * pedgerdiv / rhoedger * pedgergra * 0.25 * met13edger / dz
    else
      ardd = 0.0
    end

    # ------------------ A(i-1,j,k+2) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      aluu =
        jacinv / dx * pedgeldiv / rhoedgel * pedgelgra * 0.25 * met13edgel / dz
    else
      aluu = 0.0
    end

    # ------------------ A(i-1,j,k-2) --------------------

    if k == k1 && zboundaries == SolidWallBoundaries()
      aldd =
        -jacinv / dx * pedgeldiv / rhoedgel * pedgelgra * 0.25 * met13edgel / dz
    else
      aldd = 0.0
    end

    # ------------------ A(i,j+1,k+2) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      afuu =
        -jacinv / dy * pedgefdiv / rhoedgef * pedgefgra * 0.25 * met23edgef / dz
    else
      afuu = 0.0
    end

    # ------------------ A(i,j+1,k-2) --------------------

    if k == k1 && zboundaries == SolidWallBoundaries()
      afdd =
        jacinv / dy * pedgefdiv / rhoedgef * pedgefgra * 0.25 * met23edgef / dz
    else
      afdd = 0.0
    end

    # ------------------ A(i,j-1,k+2) --------------------

    if k == k0 && zboundaries == SolidWallBoundaries()
      abuu =
        jacinv / dy * pedgebdiv / rhoedgeb * pedgebgra * 0.25 * met23edgeb / dz
    else
      abuu = 0.0
    end

    # ------------------ A(i,j-1,k-2) --------------------

    if k == k1 && zboundaries == SolidWallBoundaries()
      abdd =
        -jacinv / dy * pedgebdiv / rhoedgeb * pedgebgra * 0.25 * met23edgeb / dz
    else
      abdd = 0.0
    end

    # Scale the tensor elements.
    ac = ac / (fcscal^2.0)
    ar = ar / fcscal / fcscal_r
    al = al / fcscal / fcscal_l
    af = af / fcscal / fcscal_f
    ab = ab / fcscal / fcscal_b
    au = au / fcscal / fcscal_u
    ad = ad / fcscal / fcscal_d
    aru = aru / fcscal / fcscal_ru
    ard = ard / fcscal / fcscal_rd
    alu = alu / fcscal / fcscal_lu
    ald = ald / fcscal / fcscal_ld
    afu = afu / fcscal / fcscal_fu
    afd = afd / fcscal / fcscal_fd
    abu = abu / fcscal / fcscal_bu
    abd = abd / fcscal / fcscal_bd
    auu = auu / fcscal / fcscal_uu
    add = add / fcscal / fcscal_dd
    aruu = aruu / fcscal / fcscal_ruu
    ardd = ardd / fcscal / fcscal_rdd
    aluu = aluu / fcscal / fcscal_luu
    aldd = aldd / fcscal / fcscal_ldd
    afuu = afuu / fcscal / fcscal_fuu
    afdd = afdd / fcscal / fcscal_fdd
    abuu = abuu / fcscal / fcscal_buu
    abdd = abdd / fcscal / fcscal_bdd

    # Determine indices for the operator.
    ia = i - i0 + 1
    ja = j - j0 + 1
    ka = k - k0 + 1

    # Set tensor elements for bicgstab.
    ac_b[ia, ja, ka] = ac
    ar_b[ia, ja, ka] = ar
    al_b[ia, ja, ka] = al
    af_b[ia, ja, ka] = af
    ab_b[ia, ja, ka] = ab
    au_b[ia, ja, ka] = au
    ad_b[ia, ja, ka] = ad
    aru_b[ia, ja, ka] = aru
    ard_b[ia, ja, ka] = ard
    alu_b[ia, ja, ka] = alu
    ald_b[ia, ja, ka] = ald
    afu_b[ia, ja, ka] = afu
    afd_b[ia, ja, ka] = afd
    abu_b[ia, ja, ka] = abu
    abd_b[ia, ja, ka] = abd
    auu_b[ia, ja, ka] = auu
    add_b[ia, ja, ka] = add
    aruu_b[ia, ja, ka] = aruu
    ardd_b[ia, ja, ka] = ardd
    aluu_b[ia, ja, ka] = aluu
    aldd_b[ia, ja, ka] = aldd
    afuu_b[ia, ja, ka] = afuu
    afdd_b[ia, ja, ka] = afdd
    abuu_b[ia, ja, ka] = abuu
    abdd_b[ia, ja, ka] = abdd

    # Store horizontal and vertical components of AC (for
    # preconditioner).
    if preconditioner
      ach_b[ia, ja, ka] = -ar - al - af - ab
      acv_b[ia, ja, ka] = -au - ad
    end
  end
  return
end

function compute_operator!(
  state::State,
  dt::AbstractFloat,
  opt::IMPL,
  facray::AbstractFloat,
)
  (; preconditioner) = state.namelists.poisson
  (; spongelayer, sponge_uv) = state.namelists.sponge
  (; zboundaries) = state.namelists.setting
  (; i0, i1, j0, j1, k0, k1) = state.domain
  (; dx, dy, dz, jac, met) = state.grid
  (; pstrattfc, rhostrattfc, bvsstrattfc) = state.atmosphere
  (; kr_sp_tfc, kr_sp_w_tfc) = state.sponge
  (;
    ac_b,
    acv_b,
    ach_b,
    al_b,
    ar_b,
    ab_b,
    af_b,
    ad_b,
    au_b,
    aru_b,
    ard_b,
    alu_b,
    ald_b,
    afu_b,
    afd_b,
    abu_b,
    abd_b,
    auu_b,
    add_b,
    aruu_b,
    ardd_b,
    aluu_b,
    aldd_b,
    afuu_b,
    afdd_b,
    abuu_b,
    abdd_b,
  ) = state.poisson.tensor
  (; rho) = state.variables.predictands

  # Compute tensor elements for TFC.
  for k in k0:k1, j in j0:j1, i in i0:i1
    # Compute scaling factors.
    fcscal = sqrt(pstrattfc[i, j, k]^2.0 / rhostrattfc[i, j, k])
    fcscal_r = sqrt(pstrattfc[i + 1, j, k]^2.0 / rhostrattfc[i + 1, j, k])
    fcscal_l = sqrt(pstrattfc[i - 1, j, k]^2.0 / rhostrattfc[i - 1, j, k])
    fcscal_f = sqrt(pstrattfc[i, j + 1, k]^2.0 / rhostrattfc[i, j + 1, k])
    fcscal_b = sqrt(pstrattfc[i, j - 1, k]^2.0 / rhostrattfc[i, j - 1, k])
    fcscal_u = sqrt(pstrattfc[i, j, k + 1]^2.0 / rhostrattfc[i, j, k + 1])
    fcscal_d = sqrt(pstrattfc[i, j, k - 1]^2.0 / rhostrattfc[i, j, k - 1])
    fcscal_ru =
      sqrt(pstrattfc[i + 1, j, k + 1]^2.0 / rhostrattfc[i + 1, j, k + 1])
    fcscal_rd =
      sqrt(pstrattfc[i + 1, j, k - 1]^2.0 / rhostrattfc[i + 1, j, k - 1])
    fcscal_lu =
      sqrt(pstrattfc[i - 1, j, k + 1]^2.0 / rhostrattfc[i - 1, j, k + 1])
    fcscal_ld =
      sqrt(pstrattfc[i - 1, j, k - 1]^2.0 / rhostrattfc[i - 1, j, k - 1])
    fcscal_fu =
      sqrt(pstrattfc[i, j + 1, k + 1]^2.0 / rhostrattfc[i, j + 1, k + 1])
    fcscal_fd =
      sqrt(pstrattfc[i, j + 1, k - 1]^2.0 / rhostrattfc[i, j + 1, k - 1])
    fcscal_bu =
      sqrt(pstrattfc[i, j - 1, k + 1]^2.0 / rhostrattfc[i, j - 1, k + 1])
    fcscal_bd =
      sqrt(pstrattfc[i, j - 1, k - 1]^2.0 / rhostrattfc[i, j - 1, k - 1])
    fcscal_uu = sqrt(pstrattfc[i, j, k + 2]^2.0 / rhostrattfc[i, j, k + 2])
    fcscal_dd = sqrt(pstrattfc[i, j, k - 2]^2.0 / rhostrattfc[i, j, k - 2])
    fcscal_ruu =
      sqrt(pstrattfc[i + 1, j, k + 2]^2.0 / rhostrattfc[i + 1, j, k + 2])
    fcscal_rdd =
      sqrt(pstrattfc[i + 1, j, k - 2]^2.0 / rhostrattfc[i + 1, j, k - 2])
    fcscal_luu =
      sqrt(pstrattfc[i - 1, j, k + 2]^2.0 / rhostrattfc[i - 1, j, k + 2])
    fcscal_ldd =
      sqrt(pstrattfc[i - 1, j, k - 2]^2.0 / rhostrattfc[i - 1, j, k - 2])
    fcscal_fuu =
      sqrt(pstrattfc[i, j + 1, k + 2]^2.0 / rhostrattfc[i, j + 1, k + 2])
    fcscal_fdd =
      sqrt(pstrattfc[i, j + 1, k - 2]^2.0 / rhostrattfc[i, j + 1, k - 2])
    fcscal_buu =
      sqrt(pstrattfc[i, j - 1, k + 2]^2.0 / rhostrattfc[i, j - 1, k + 2])
    fcscal_bdd =
      sqrt(pstrattfc[i, j - 1, k - 2]^2.0 / rhostrattfc[i, j - 1, k - 2])

    # Compute inverse Jacobian.
    jacinv = 1.0 / jac[i, j, k]

    # Compute P coefficients (divergence).
    pedgerdiv =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
      )
    pedgeldiv =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i - 1, j, k] * pstrattfc[i - 1, j, k]
      )
    pedgefdiv =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
      )
    pedgebdiv =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j - 1, k] * pstrattfc[i, j - 1, k]
      )
    pedgeudiv =
      jac[i, j, k] *
      jac[i, j, k + 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1])
    pedgeddiv =
      jac[i, j, k] *
      jac[i, j, k - 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k - 1]) /
      (jac[i, j, k] + jac[i, j, k - 1])

    # Compute P coefficients (pressure gradient).
    pedgergra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
    pedgelgra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i - 1, j, k])
    pedgefgra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
    pedgebgra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j - 1, k])
    pedgeugra =
      (
        jac[i, j, k + 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    pedgedgra =
      (
        jac[i, j, k - 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k - 1]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    puedgergra = 0.5 * (pstrattfc[i, j, k + 1] + pstrattfc[i + 1, j, k + 1])
    puedgelgra = 0.5 * (pstrattfc[i, j, k + 1] + pstrattfc[i - 1, j, k + 1])
    puedgefgra = 0.5 * (pstrattfc[i, j, k + 1] + pstrattfc[i, j + 1, k + 1])
    puedgebgra = 0.5 * (pstrattfc[i, j, k + 1] + pstrattfc[i, j - 1, k + 1])
    pdedgergra = 0.5 * (pstrattfc[i, j, k - 1] + pstrattfc[i + 1, j, k - 1])
    pdedgelgra = 0.5 * (pstrattfc[i, j, k - 1] + pstrattfc[i - 1, j, k - 1])
    pdedgefgra = 0.5 * (pstrattfc[i, j, k - 1] + pstrattfc[i, j + 1, k - 1])
    pdedgebgra = 0.5 * (pstrattfc[i, j, k - 1] + pstrattfc[i, j - 1, k - 1])

    # Compute density coefficients.
    rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
    rhostratedgel = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i - 1, j, k])
    rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
    rhostratedgeb = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j - 1, k])
    rhostratedgeu =
      (
        jac[i, j, k + 1] * rhostrattfc[i, j, k] +
        jac[i, j, k] * rhostrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    rhostratedged =
      (
        jac[i, j, k - 1] * rhostrattfc[i, j, k] +
        jac[i, j, k] * rhostrattfc[i, j, k - 1]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    rhoedger = 0.5 * (rho[i, j, k] + rho[i + 1, j, k]) + rhostratedger
    rhoedgel = 0.5 * (rho[i, j, k] + rho[i - 1, j, k]) + rhostratedgel
    rhoedgef = 0.5 * (rho[i, j, k] + rho[i, j + 1, k]) + rhostratedgef
    rhoedgeb = 0.5 * (rho[i, j, k] + rho[i, j - 1, k]) + rhostratedgeb
    rhoedgeu =
      (jac[i, j, k + 1] * rho[i, j, k] + jac[i, j, k] * rho[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1]) + rhostratedgeu
    rhoedged =
      (jac[i, j, k - 1] * rho[i, j, k] + jac[i, j, k] * rho[i, j, k - 1]) /
      (jac[i, j, k] + jac[i, j, k - 1]) + rhostratedged

    rhouedger =
      0.5 * (
        rho[i, j, k + 1] +
        rho[i + 1, j, k + 1] +
        rhostrattfc[i, j, k + 1] +
        rhostrattfc[i + 1, j, k + 1]
      )
    rhouedgel =
      0.5 * (
        rho[i, j, k + 1] +
        rho[i - 1, j, k + 1] +
        rhostrattfc[i, j, k + 1] +
        rhostrattfc[i - 1, j, k + 1]
      )
    rhouedgef =
      0.5 * (
        rho[i, j, k + 1] +
        rho[i, j + 1, k + 1] +
        rhostrattfc[i, j, k + 1] +
        rhostrattfc[i, j + 1, k + 1]
      )
    rhouedgeb =
      0.5 * (
        rho[i, j, k + 1] +
        rho[i, j - 1, k + 1] +
        rhostrattfc[i, j, k + 1] +
        rhostrattfc[i, j - 1, k + 1]
      )
    rhodedger =
      0.5 * (
        rho[i, j, k - 1] +
        rho[i + 1, j, k - 1] +
        rhostrattfc[i, j, k - 1] +
        rhostrattfc[i + 1, j, k - 1]
      )
    rhodedgel =
      0.5 * (
        rho[i, j, k - 1] +
        rho[i - 1, j, k - 1] +
        rhostrattfc[i, j, k - 1] +
        rhostrattfc[i - 1, j, k - 1]
      )
    rhodedgef =
      0.5 * (
        rho[i, j, k - 1] +
        rho[i, j + 1, k - 1] +
        rhostrattfc[i, j, k - 1] +
        rhostrattfc[i, j + 1, k - 1]
      )
    rhodedgeb =
      0.5 * (
        rho[i, j, k - 1] +
        rho[i, j - 1, k - 1] +
        rhostrattfc[i, j, k - 1] +
        rhostrattfc[i, j - 1, k - 1]
      )

    # Compute squared buoyancy frequency at edges.
    bvsstratedgeu =
      (
        jac[i, j, k + 1] * bvsstrattfc[i, j, k] +
        jac[i, j, k] * bvsstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    bvsstratedged =
      (
        jac[i, j, k - 1] * bvsstrattfc[i, j, k] +
        jac[i, j, k] * bvsstrattfc[i, j, k - 1]
      ) / (jac[i, j, k] + jac[i, j, k - 1])

    # Interpolate metric-tensor elements.
    met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])
    met13edgel = 0.5 * (met[i, j, k, 1, 3] + met[i - 1, j, k, 1, 3])
    met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])
    met23edgeb = 0.5 * (met[i, j, k, 2, 3] + met[i, j - 1, k, 2, 3])
    met13edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 1, 3] +
        jac[i, j, k] * met[i, j, k + 1, 1, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met23edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 2, 3] +
        jac[i, j, k] * met[i, j, k + 1, 2, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met33edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 3, 3] +
        jac[i, j, k] * met[i, j, k + 1, 3, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met13edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 1, 3] +
        jac[i, j, k] * met[i, j, k - 1, 1, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    met23edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 2, 3] +
        jac[i, j, k] * met[i, j, k - 1, 2, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    met33edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 3, 3] +
        jac[i, j, k] * met[i, j, k - 1, 3, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    met13uedger = 0.5 * (met[i, j, k + 1, 1, 3] + met[i + 1, j, k + 1, 1, 3])
    met13uedgel = 0.5 * (met[i, j, k + 1, 1, 3] + met[i - 1, j, k + 1, 1, 3])
    met23uedgef = 0.5 * (met[i, j, k + 1, 2, 3] + met[i, j + 1, k + 1, 2, 3])
    met23uedgeb = 0.5 * (met[i, j, k + 1, 2, 3] + met[i, j - 1, k + 1, 2, 3])
    met13dedger = 0.5 * (met[i, j, k - 1, 1, 3] + met[i + 1, j, k - 1, 1, 3])
    met13dedgel = 0.5 * (met[i, j, k - 1, 1, 3] + met[i - 1, j, k - 1, 1, 3])
    met23dedgef = 0.5 * (met[i, j, k - 1, 2, 3] + met[i, j + 1, k - 1, 2, 3])
    met23dedgeb = 0.5 * (met[i, j, k - 1, 2, 3] + met[i, j - 1, k - 1, 2, 3])

    # Compute Rayleigh damping terms.
    facedger = 1.0
    facedgel = 1.0
    facedgef = 1.0
    facedgeb = 1.0
    facuedger = 1.0
    facuedgel = 1.0
    facuedgef = 1.0
    facuedgeb = 1.0
    facdedger = 1.0
    facdedgel = 1.0
    facdedgef = 1.0
    facdedgeb = 1.0
    facedgeu = 1.0
    facedged = 1.0
    if spongelayer
      if sponge_uv
        facedger =
          facedger +
          dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i + 1, j, k]) * facray
        facedgel =
          facedgel +
          dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i - 1, j, k]) * facray
        facedgef =
          facedgef +
          dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j + 1, k]) * facray
        facedgeb =
          facedgeb +
          dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j - 1, k]) * facray
        facuedger =
          facuedger +
          dt *
          0.5 *
          (kr_sp_tfc[i, j, k + 1] + kr_sp_tfc[i + 1, j, k + 1]) *
          facray
        facuedgel =
          facuedgel +
          dt *
          0.5 *
          (kr_sp_tfc[i, j, k + 1] + kr_sp_tfc[i - 1, j, k + 1]) *
          facray
        facuedgef =
          facuedgef +
          dt *
          0.5 *
          (kr_sp_tfc[i, j, k + 1] + kr_sp_tfc[i, j + 1, k + 1]) *
          facray
        facuedgeb =
          facuedgeb +
          dt *
          0.5 *
          (kr_sp_tfc[i, j, k + 1] + kr_sp_tfc[i, j - 1, k + 1]) *
          facray
        facdedger =
          facdedger +
          dt *
          0.5 *
          (kr_sp_tfc[i, j, k - 1] + kr_sp_tfc[i + 1, j, k - 1]) *
          facray
        facdedgel =
          facdedgel +
          dt *
          0.5 *
          (kr_sp_tfc[i, j, k - 1] + kr_sp_tfc[i - 1, j, k - 1]) *
          facray
        facdedgef =
          facdedgef +
          dt *
          0.5 *
          (kr_sp_tfc[i, j, k - 1] + kr_sp_tfc[i, j + 1, k - 1]) *
          facray
        facdedgeb =
          facdedgeb +
          dt *
          0.5 *
          (kr_sp_tfc[i, j, k - 1] + kr_sp_tfc[i, j - 1, k - 1]) *
          facray
      end
      facedgeu =
        facedgeu +
        dt * (
          jac[i, j, k + 1] * kr_sp_w_tfc[i, j, k] +
          jac[i, j, k] * kr_sp_w_tfc[i, j, k + 1]
        ) / (jac[i, j, k] + jac[i, j, k + 1]) * facray
      facedged =
        facedged +
        dt * (
          jac[i, j, k - 1] * kr_sp_w_tfc[i, j, k] +
          jac[i, j, k] * kr_sp_w_tfc[i, j, k - 1]
        ) / (jac[i, j, k] + jac[i, j, k - 1]) * facray
    end

    # Compute implicit coefficients.
    imphoredger = 1.0 / (facedger^2.0)
    imphoredgel = 1.0 / (facedgel^2.0)
    imphoredgef = 1.0 / (facedgef^2.0)
    imphoredgeb = 1.0 / (facedgeb^2.0)
    imphoruedger = 1.0 / (facuedger^2.0)
    imphoruedgel = 1.0 / (facuedgel^2.0)
    imphoruedgef = 1.0 / (facuedgef^2.0)
    imphoruedgeb = 1.0 / (facuedgeb^2.0)
    imphordedger = 1.0 / (facdedger^2.0)
    imphordedgel = 1.0 / (facdedgel^2.0)
    imphordedgef = 1.0 / (facdedgef^2.0)
    imphordedgeb = 1.0 / (facdedgeb^2.0)
    impveredgeu =
      1.0 / (facedgeu + rhostratedgeu / rhoedgeu * bvsstratedgeu * dt^2.0)
    impveredged =
      1.0 / (facedged + rhostratedged / rhoedged * bvsstratedged * dt^2.0)

    # Compute gradient coefficients

    # G(i + 1 / 2)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gedger =
        jacinv / dx * pedgerdiv * imphoredger * facedger / rhoedger +
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k, 1, 3] *
        jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoredger *
        facedger / rhoedger
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      gedger =
        jacinv / dx * pedgerdiv * imphoredger * facedger / rhoedger -
        jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k, 1, 3] *
        jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphoredger *
        facedger / rhoedger
    else
      gedger =
        jacinv / dx * pedgerdiv * imphoredger * facedger / rhoedger +
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k, 1, 3] *
        jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoredger *
        facedger / rhoedger -
        jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k, 1, 3] *
        jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphoredger *
        facedger / rhoedger
    end

    # G(i - 1 / 2)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gedgel =
        -jacinv / dx * pedgeldiv * imphoredgel * facedgel / rhoedgel +
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k, 1, 3] *
        jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoredgel *
        facedgel / rhoedgel
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      gedgel =
        -jacinv / dx * pedgeldiv * imphoredgel * facedgel / rhoedgel -
        jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k, 1, 3] *
        jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphoredgel *
        facedgel / rhoedgel
    else
      gedgel =
        -jacinv / dx * pedgeldiv * imphoredgel * facedgel / rhoedgel +
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k, 1, 3] *
        jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoredgel *
        facedgel / rhoedgel -
        jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k, 1, 3] *
        jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphoredgel *
        facedgel / rhoedgel
    end

    # G(j + 1 / 2)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gedgef =
        jacinv / dy * pedgefdiv * imphoredgef * facedgef / rhoedgef +
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k, 2, 3] *
        jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoredgef *
        facedgef / rhoedgef
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      gedgef =
        jacinv / dy * pedgefdiv * imphoredgef * facedgef / rhoedgef -
        jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k, 2, 3] *
        jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphoredgef *
        facedgef / rhoedgef
    else
      gedgef =
        jacinv / dy * pedgefdiv * imphoredgef * facedgef / rhoedgef +
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k, 2, 3] *
        jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoredgef *
        facedgef / rhoedgef -
        jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k, 2, 3] *
        jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphoredgef *
        facedgef / rhoedgef
    end

    # G(j - 1 / 2)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gedgeb =
        -jacinv / dy * pedgebdiv * imphoredgeb * facedgeb / rhoedgeb +
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k, 2, 3] *
        jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoredgeb *
        facedgeb / rhoedgeb
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      gedgeb =
        -jacinv / dy * pedgebdiv * imphoredgeb * facedgeb / rhoedgeb -
        jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k, 2, 3] *
        jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphoredgeb *
        facedgeb / rhoedgeb
    else
      gedgeb =
        -jacinv / dy * pedgebdiv * imphoredgeb * facedgeb / rhoedgeb +
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k, 2, 3] *
        jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoredgeb *
        facedgeb / rhoedgeb -
        jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k, 2, 3] *
        jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphoredgeb *
        facedgeb / rhoedgeb
    end

    # G(k + 1 / 2)
    if k == k1 && zboundaries == SolidWallBoundaries()
      gedgeu = 0.0
    else
      gedgeu = jacinv / dz * pedgeudiv * impveredgeu / rhoedgeu
    end

    # G(k - 1 / 2)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gedged = 0.0
    else
      gedged = -jacinv / dz * pedgeddiv * impveredged / rhoedged
    end

    # G(i + 1 / 2, k + 1)
    if k == k1 && zboundaries == SolidWallBoundaries()
      guedger = 0.0
    else
      guedger =
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k + 1, 1, 3] *
        jac[i, j, k] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoruedger *
        facuedger / rhouedger
    end

    # G(i - 1 / 2, k + 1)
    if k == k1 && zboundaries == SolidWallBoundaries()
      guedgel = 0.0
    else
      guedgel =
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k + 1, 1, 3] *
        jac[i, j, k] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoruedgel *
        facuedgel / rhouedgel
    end

    # G(j + 1 / 2, k + 1)
    if k == k1 && zboundaries == SolidWallBoundaries()
      guedgef = 0.0
    else
      guedgef =
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k + 1, 2, 3] *
        jac[i, j, k] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoruedgef *
        facuedgef / rhouedgef
    end

    # G(j - 1 / 2, k + 1)
    if k == k1 && zboundaries == SolidWallBoundaries()
      guedgeb = 0.0
    else
      guedgeb =
        jacinv / dz * pedgeudiv * impveredgeu * rhostratedgeu / rhoedgeu *
        bvsstratedgeu *
        dt^2.0 *
        0.5 *
        met[i, j, k + 1, 2, 3] *
        jac[i, j, k] / (jac[i, j, k] + jac[i, j, k + 1]) *
        imphoruedgeb *
        facuedgeb / rhouedgeb
    end

    # G(i + 1 / 2, k - 1)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gdedger = 0.0
    else
      gdedger =
        -jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k - 1, 1, 3] *
        jac[i, j, k] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphordedger *
        facdedger / rhodedger
    end

    # G(i - 1 / 2, k - 1)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gdedgel = 0.0
    else
      gdedgel =
        -jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k - 1, 1, 3] *
        jac[i, j, k] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphordedgel *
        facdedgel / rhodedgel
    end

    # G(j + 1 / 2, k - 1)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gdedgef = 0.0
    else
      gdedgef =
        -jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k - 1, 2, 3] *
        jac[i, j, k] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphordedgef *
        facdedgef / rhodedgef
    end

    # G(j - 1 / 2, k - 1)
    if k == k0 && zboundaries == SolidWallBoundaries()
      gdedgeb = 0.0
    else
      gdedgeb =
        -jacinv / dz * pedgeddiv * impveredged * rhostratedged / rhoedged *
        bvsstratedged *
        dt^2.0 *
        0.5 *
        met[i, j, k - 1, 2, 3] *
        jac[i, j, k] / (jac[i, j, k] + jac[i, j, k - 1]) *
        imphordedgeb *
        facdedgeb / rhodedgeb
    end

    # Compute tensor elements

    # ------------------- A(i,j,k) --------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      ac =
        -gedger * pedgergra * (1.0 / dx + 0.75 * met13edger / dz) +
        gedgel * pedgelgra * (1.0 / dx - 0.75 * met13edgel / dz) -
        gedgef * pedgefgra * (1.0 / dy + 0.75 * met23edgef / dz) +
        gedgeb * pedgebgra * (1.0 / dy - 0.75 * met23edgeb / dz) -
        gedgeu * pedgeugra * met33edgeu / dz -
        guedger * puedgergra * 0.25 * met13uedger / dz -
        guedgel * puedgelgra * 0.25 * met13uedgel / dz -
        guedgef * puedgefgra * 0.25 * met23uedgef / dz -
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      ac =
        -gedger * pedgergra / dx + gedgel * pedgelgra / dx -
        gedgef * pedgefgra / dy + gedgeb * pedgebgra / dy -
        gedgeu * pedgeugra * met33edgeu / dz +
        gedged * pedgedgra * met33edged / dz -
        guedger * puedgergra * 0.25 * met13uedger / dz -
        guedgel * puedgelgra * 0.25 * met13uedgel / dz -
        guedgef * puedgefgra * 0.25 * met23uedgef / dz -
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz +
        gdedger * pdedgergra * met13dedger / dz +
        gdedgel * pdedgelgra * met13dedgel / dz +
        gdedgef * pdedgefgra * met23dedgef / dz +
        gdedgeb * pdedgebgra * met23dedgeb / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      ac =
        -gedger * pedgergra / dx + gedgel * pedgelgra / dx -
        gedgef * pedgefgra / dy + gedgeb * pedgebgra / dy -
        gedgeu * pedgeugra * met33edgeu / dz +
        gedged * pedgedgra * met33edged / dz -
        guedger * puedgergra * met13uedger / dz -
        guedgel * puedgelgra * met13uedgel / dz -
        guedgef * puedgefgra * met23uedgef / dz -
        guedgeb * puedgebgra * met23uedgeb / dz +
        gdedger * pdedgergra * 0.25 * met13dedger / dz +
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz +
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz +
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ac =
        -gedger * pedgergra * (1.0 / dx - 0.75 * met13edger / dz) +
        gedgel * pedgelgra * (1.0 / dx + 0.75 * met13edgel / dz) -
        gedgef * pedgefgra * (1.0 / dy - 0.75 * met23edgef / dz) +
        gedgeb * pedgebgra * (1.0 / dy + 0.75 * met23edgeb / dz) +
        gedged * pedgedgra * met33edged / dz +
        gdedger * pdedgergra * 0.25 * met13dedger / dz +
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz +
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz +
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    else
      ac =
        -gedger * pedgergra / dx + gedgel * pedgelgra / dx -
        gedgef * pedgefgra / dy + gedgeb * pedgebgra / dy -
        gedgeu * pedgeugra * met33edgeu / dz +
        gedged * pedgedgra * met33edged / dz -
        guedger * puedgergra * 0.25 * met13uedger / dz -
        guedgel * puedgelgra * 0.25 * met13uedgel / dz -
        guedgef * puedgefgra * 0.25 * met23uedgef / dz -
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz +
        gdedger * pdedgergra * 0.25 * met13dedger / dz +
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz +
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz +
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    end

    # ------------------ A(i+1,j,k) -------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      ar =
        gedger * pedgergra * (1.0 / dx - 0.75 * met13edger / dz) +
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k + 1] /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) -
        guedger * puedgergra * 0.25 * met13uedger / dz
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      ar =
        gedger * pedgergra / dx +
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k + 1] /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k - 1] /
        (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) -
        guedger * puedgergra * 0.25 * met13uedger / dz +
        gdedger * pdedgergra * met13dedger / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      ar =
        gedger * pedgergra / dx +
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k + 1] /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k - 1] /
        (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) -
        guedger * puedgergra * met13uedger / dz +
        gdedger * pdedgergra * 0.25 * met13dedger / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ar =
        gedger * pedgergra * (1.0 / dx + 0.75 * met13edger / dz) +
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k - 1] /
        (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
        gdedger * pdedgergra * 0.25 * met13dedger / dz
    else
      ar =
        gedger * pedgergra / dx +
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k + 1] /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k - 1] /
        (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) -
        guedger * puedgergra * 0.25 * met13uedger / dz +
        gdedger * pdedgergra * 0.25 * met13dedger / dz
    end

    # ------------------ A(i-1,j,k) -------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      al =
        -gedgel * pedgelgra * (1.0 / dx + 0.75 * met13edgel / dz) -
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k + 1] /
        (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
        guedgel * puedgelgra * 0.25 * met13uedgel / dz
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      al =
        -gedgel * pedgelgra / dx -
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k + 1] /
        (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k - 1] /
        (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
        guedgel * puedgelgra * 0.25 * met13uedgel / dz +
        gdedgel * pdedgelgra * met13dedgel / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      al =
        -gedgel * pedgelgra / dx -
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k + 1] /
        (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k - 1] /
        (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
        guedgel * puedgelgra * met13uedgel / dz +
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      al =
        -gedgel * pedgelgra * (1.0 / dx - 0.75 * met13edgel / dz) -
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k - 1] /
        (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) +
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
    else
      al =
        -gedgel * pedgelgra / dx -
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k + 1] /
        (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k - 1] /
        (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
        guedgel * puedgelgra * 0.25 * met13uedgel / dz +
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
    end

    # ------------------ A(i,j+1,k) -------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      af =
        gedgef * pedgefgra * (1.0 / dy - 0.75 * met23edgef / dz) +
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k + 1] /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) -
        guedgef * puedgefgra * 0.25 * met23uedgef / dz
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      af =
        gedgef * pedgefgra / dy +
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k + 1] /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k - 1] /
        (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) -
        guedgef * puedgefgra * 0.25 * met23uedgef / dz +
        gdedgef * pdedgefgra * met23dedgef / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      af =
        gedgef * pedgefgra / dy +
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k + 1] /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k - 1] /
        (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) -
        guedgef * puedgefgra * met23uedgef / dz +
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      af =
        gedgef * pedgefgra * (1.0 / dy + 0.75 * met23edgef / dz) +
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k - 1] /
        (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
    else
      af =
        gedgef * pedgefgra / dy +
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k + 1] /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k - 1] /
        (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) -
        guedgef * puedgefgra * 0.25 * met23uedgef / dz +
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
    end

    # ------------------ A(i,j-1,k) -------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      ab =
        -gedgeb * pedgebgra * (1.0 / dy + 0.75 * met23edgeb / dz) -
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k + 1] /
        (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      ab =
        -gedgeb * pedgebgra / dy -
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k + 1] /
        (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k - 1] /
        (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz +
        gdedgeb * pdedgebgra * met23dedgeb / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      ab =
        -gedgeb * pedgebgra / dy -
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k + 1] /
        (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k - 1] /
        (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
        guedgeb * puedgebgra * met23uedgeb / dz +
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ab =
        -gedgeb * pedgebgra * (1.0 / dy - 0.75 * met23edgeb / dz) -
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k - 1] /
        (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) +
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    else
      ab =
        -gedgeb * pedgebgra / dy -
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k + 1] /
        (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k - 1] /
        (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz +
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    end

    # ------------------ A(i,j,k+1) -------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      au =
        gedger * pedgergra * met13edger / dz +
        gedgel * pedgelgra * met13edgel / dz +
        gedgef * pedgefgra * met23edgef / dz +
        gedgeb * pedgebgra * met23edgeb / dz +
        gedgeu * pedgeugra * met33edgeu / dz - guedger * puedgergra / dx +
        guedgel * puedgelgra / dx - guedgef * puedgefgra / dy +
        guedgeb * puedgebgra / dy
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      au =
        gedger * pedgergra * 0.25 * met13edger / dz +
        gedgel * pedgelgra * 0.25 * met13edgel / dz +
        gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
        gedgeu * pedgeugra * met33edgeu / dz - guedger * puedgergra / dx +
        guedgel * puedgelgra / dx - guedgef * puedgefgra / dy +
        guedgeb * puedgebgra / dy -
        gdedger * pdedgergra * 0.25 * met13dedger / dz -
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz -
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz -
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      au =
        gedger * pedgergra * 0.25 * met13edger / dz +
        gedgel * pedgelgra * 0.25 * met13edgel / dz +
        gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
        gedgeu * pedgeugra * met33edgeu / dz -
        guedger * puedgergra * (1.0 / dx - 0.75 * met13uedger / dz) +
        guedgel * puedgelgra * (1.0 / dx + 0.75 * met13uedgel / dz) -
        guedgef * puedgefgra * (1.0 / dy - 0.75 * met23uedgef / dz) +
        guedgeb * puedgebgra * (1.0 / dy + 0.75 * met23uedgeb / dz)
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      au = 0.0
    else
      au =
        gedger * pedgergra * 0.25 * met13edger / dz +
        gedgel * pedgelgra * 0.25 * met13edgel / dz +
        gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
        gedgeu * pedgeugra * met33edgeu / dz - guedger * puedgergra / dx +
        guedgel * puedgelgra / dx - guedgef * puedgefgra / dy +
        guedgeb * puedgebgra / dy
    end

    # ------------------ A(i,j,k-1) -------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      ad = 0.0
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      ad =
        -gedger * pedgergra * 0.25 * met13edger / dz -
        gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedgef * pedgefgra * 0.25 * met23edgef / dz -
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedged * pedgedgra * met33edged / dz -
        gdedger * pdedgergra * (1.0 / dx + 0.75 * met13dedger / dz) +
        gdedgel * pdedgelgra * (1.0 / dx - 0.75 * met13dedgel / dz) -
        gdedgef * pdedgefgra * (1.0 / dy + 0.75 * met23dedgef / dz) +
        gdedgeb * pdedgebgra * (1.0 / dy - 0.75 * met23dedgeb / dz)
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      ad =
        -gedger * pedgergra * 0.25 * met13edger / dz -
        gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedgef * pedgefgra * 0.25 * met23edgef / dz -
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedged * pedgedgra * met33edged / dz - gdedger * pdedgergra / dx +
        gdedgel * pdedgelgra / dx - gdedgef * pdedgefgra / dy +
        gdedgeb * pdedgebgra / dy +
        guedger * puedgergra * 0.25 * met13uedger / dz +
        guedgel * puedgelgra * 0.25 * met13uedgel / dz +
        guedgef * puedgefgra * 0.25 * met23uedgef / dz +
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ad =
        -gedger * pedgergra * met13edger / dz -
        gedgel * pedgelgra * met13edgel / dz -
        gedgef * pedgefgra * met23edgef / dz -
        gedgeb * pedgebgra * met23edgeb / dz -
        gedged * pedgedgra * met33edged / dz - gdedger * pdedgergra / dx +
        gdedgel * pdedgelgra / dx - gdedgef * pdedgefgra / dy +
        gdedgeb * pdedgebgra / dy
    else
      ad =
        -gedger * pedgergra * 0.25 * met13edger / dz -
        gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedgef * pedgefgra * 0.25 * met23edgef / dz -
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedged * pedgedgra * met33edged / dz - gdedger * pdedgergra / dx +
        gdedgel * pdedgelgra / dx - gdedgef * pdedgefgra / dy +
        gdedgeb * pdedgebgra / dy
    end

    # ----------------- A(i+1,j,k+1) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      aru =
        gedger * pedgergra * met13edger / dz +
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k] /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
        guedger * puedgergra / dx
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      aru =
        gedger * pedgergra * 0.25 * met13edger / dz +
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k] /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
        guedger * puedgergra / dx -
        gdedger * pdedgergra * 0.25 * met13dedger / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      aru =
        gedger * pedgergra * 0.25 * met13edger / dz +
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k] /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
        guedger * puedgergra * (1.0 / dx + 0.75 * met13uedger / dz)
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      aru = 0.0
    else
      aru =
        gedger * pedgergra * 0.25 * met13edger / dz +
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k] /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
        guedger * puedgergra / dx
    end

    # ----------------- A(i+1,j,k-1) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      ard = 0.0
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      ard =
        -gedger * pedgergra * 0.25 * met13edger / dz +
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k] /
        (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
        gdedger * pdedgergra * (1.0 / dx - 0.75 * met13dedger / dz)
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      ard =
        -gedger * pedgergra * 0.25 * met13edger / dz +
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k] /
        (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
        gdedger * pdedgergra / dx +
        guedger * puedgergra * 0.25 * met13uedger / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ard =
        -gedger * pedgergra * met13edger / dz +
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k] /
        (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
        gdedger * pdedgergra / dx
    else
      ard =
        -gedger * pedgergra * 0.25 * met13edger / dz +
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k] /
        (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
        gdedger * pdedgergra / dx
    end

    # ----------------- A(i-1,j,k+1) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      alu =
        gedgel * pedgelgra * met13edgel / dz -
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k] /
        (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) - guedgel * puedgelgra / dx
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      alu =
        gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k] /
        (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) - guedgel * puedgelgra / dx -
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      alu =
        gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k] /
        (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
        guedgel * puedgelgra * (1.0 / dx - 0.75 * met13uedgel / dz)
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      alu = 0.0
    else
      alu =
        gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k] /
        (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) - guedgel * puedgelgra / dx
    end

    # ----------------- A(i-1,j,k-1) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      ald = 0.0
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      ald =
        -gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k] /
        (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
        gdedgel * pdedgelgra * (1.0 / dx + 0.75 * met13dedgel / dz)
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      ald =
        -gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k] /
        (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) - gdedgel * pdedgelgra / dx +
        guedgel * puedgelgra * 0.25 * met13uedgel / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ald =
        -gedgel * pedgelgra * met13edgel / dz -
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k] /
        (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) - gdedgel * pdedgelgra / dx
    else
      ald =
        -gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k] /
        (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) - gdedgel * pdedgelgra / dx
    end

    # ----------------- A(i,j+1,k+1) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      afu =
        gedgef * pedgefgra * met23edgef / dz +
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k] /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
        guedgef * puedgefgra / dy
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      afu =
        gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k] /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
        guedgef * puedgefgra / dy -
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      afu =
        gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k] /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
        guedgef * puedgefgra * (1.0 / dy + 0.75 * met23uedgef / dz)
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      afu = 0.0
    else
      afu =
        gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k] /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
        guedgef * puedgefgra / dy
    end

    # ----------------- A(i,j+1,k-1) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      afd = 0.0
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      afd =
        -gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k] /
        (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
        gdedgef * pdedgefgra * (1.0 / dy - 0.75 * met23dedgef / dz)
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      afd =
        -gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k] /
        (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
        gdedgef * pdedgefgra / dy +
        guedgef * puedgefgra * 0.25 * met23uedgef / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      afd =
        -gedgef * pedgefgra * met23edgef / dz +
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k] /
        (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
        gdedgef * pdedgefgra / dy
    else
      afd =
        -gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k] /
        (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
        gdedgef * pdedgefgra / dy
    end

    # ----------------- A(i,j-1,k+1) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      abu =
        gedgeb * pedgebgra * met23edgeb / dz -
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k] /
        (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) - guedgeb * puedgebgra / dy
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      abu =
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k] /
        (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) - guedgeb * puedgebgra / dy -
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      abu =
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k] /
        (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
        guedgeb * puedgebgra * (1.0 / dy - 0.75 * met23uedgeb / dz)
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      abu = 0.0
    else
      abu =
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k] /
        (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) - guedgeb * puedgebgra / dy
    end

    # ----------------- A(i,j-1,k-1) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      abd = 0.0
    elseif k == k0 + 1 && zboundaries == SolidWallBoundaries()
      abd =
        -gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k] /
        (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
        gdedgeb * pdedgebgra * (1.0 / dy + 0.75 * met23dedgeb / dz)
    elseif k == k1 - 1 && zboundaries == SolidWallBoundaries()
      abd =
        -gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k] /
        (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) - gdedgeb * pdedgebgra / dy +
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      abd =
        -gedgeb * pedgebgra * met23edgeb / dz -
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k] /
        (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) - gdedgeb * pdedgebgra / dy
    else
      abd =
        -gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k] /
        (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) - gdedgeb * pdedgebgra / dy
    end

    # ------------------ A(i,j,k+2) -------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      auu =
        -gedger * pedgergra * 0.25 * met13edger / dz -
        gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gedgef * pedgefgra * 0.25 * met23edgef / dz -
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
        guedger * puedgergra * 0.25 * met13uedger / dz +
        guedgel * puedgelgra * 0.25 * met13uedgel / dz +
        guedgef * puedgefgra * 0.25 * met23uedgef / dz +
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
    elseif (k == k1 - 1 || k == k1) && zboundaries == SolidWallBoundaries()
      auu = 0.0
    else
      auu =
        guedger * puedgergra * 0.25 * met13uedger / dz +
        guedgel * puedgelgra * 0.25 * met13uedgel / dz +
        guedgef * puedgefgra * 0.25 * met23uedgef / dz +
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
    end

    # ------------------ A(i,j,k-2) -------------------#

    if (k == k0 || k == k0 + 1) && zboundaries == SolidWallBoundaries()
      add = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      add =
        gedger * pedgergra * 0.25 * met13edger / dz +
        gedgel * pedgelgra * 0.25 * met13edgel / dz +
        gedgef * pedgefgra * 0.25 * met23edgef / dz +
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gdedger * pdedgergra * 0.25 * met13dedger / dz -
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz -
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz -
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    else
      add =
        -gdedger * pdedgergra * 0.25 * met13dedger / dz -
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz -
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz -
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    end

    # ----------------- A(i+1,j,k+2) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      aruu =
        -gedger * pedgergra * 0.25 * met13edger / dz +
        guedger * puedgergra * 0.25 * met13uedger / dz
    elseif (k == k1 - 1 || k == k1) && zboundaries == SolidWallBoundaries()
      aruu = 0.0
    else
      aruu = guedger * puedgergra * 0.25 * met13uedger / dz
    end

    # ----------------- A(i+1,j,k-2) ------------------#

    if (k == k0 || k == k0 + 1) && zboundaries == SolidWallBoundaries()
      ardd = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      ardd =
        gedger * pedgergra * 0.25 * met13edger / dz -
        gdedger * pdedgergra * 0.25 * met13dedger / dz
    else
      ardd = -gdedger * pdedgergra * 0.25 * met13dedger / dz
    end

    # ----------------- A(i-1,j,k+2) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      aluu =
        -gedgel * pedgelgra * 0.25 * met13edgel / dz +
        guedgel * puedgelgra * 0.25 * met13uedgel / dz
    elseif (k == k1 - 1 || k == k1) && zboundaries == SolidWallBoundaries()
      aluu = 0.0
    else
      aluu = guedgel * puedgelgra * 0.25 * met13uedgel / dz
    end

    # ----------------- A(i-1,j,k-2) ------------------#

    if (k == k0 || k == k0 + 1) && zboundaries == SolidWallBoundaries()
      aldd = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      aldd =
        gedgel * pedgelgra * 0.25 * met13edgel / dz -
        gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
    else
      aldd = -gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
    end

    # ----------------- A(i,j+1,k+2) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      afuu =
        -gedgef * pedgefgra * 0.25 * met23edgef / dz +
        guedgef * puedgefgra * 0.25 * met23uedgef / dz
    elseif (k == k1 - 1 || k == k1) && zboundaries == SolidWallBoundaries()
      afuu = 0.0
    else
      afuu = guedgef * puedgefgra * 0.25 * met23uedgef / dz
    end

    # ----------------- A(i,j+1,k-2) ------------------#

    if (k == k0 || k == k0 + 1) && zboundaries == SolidWallBoundaries()
      afdd = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      afdd =
        gedgef * pedgefgra * 0.25 * met23edgef / dz -
        gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
    else
      afdd = -gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
    end

    # ----------------- A(i,j-1,k+2) ------------------#

    if k == k0 && zboundaries == SolidWallBoundaries()
      abuu =
        -gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
        guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
    elseif (k == k1 - 1 || k == k1) && zboundaries == SolidWallBoundaries()
      abuu = 0.0
    else
      abuu = guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
    end

    # ----------------- A(i,j-1,k-2) ------------------#

    if (k == k0 || k == k0 + 1) && zboundaries == SolidWallBoundaries()
      abdd = 0.0
    elseif k == k1 && zboundaries == SolidWallBoundaries()
      abdd =
        gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
        gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    else
      abdd = -gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
    end

    # Scale the tensor elements.
    ac = ac / (fcscal^2.0)
    ar = ar / fcscal / fcscal_r
    al = al / fcscal / fcscal_l
    af = af / fcscal / fcscal_f
    ab = ab / fcscal / fcscal_b
    au = au / fcscal / fcscal_u
    ad = ad / fcscal / fcscal_d
    aru = aru / fcscal / fcscal_ru
    ard = ard / fcscal / fcscal_rd
    alu = alu / fcscal / fcscal_lu
    ald = ald / fcscal / fcscal_ld
    afu = afu / fcscal / fcscal_fu
    afd = afd / fcscal / fcscal_fd
    abu = abu / fcscal / fcscal_bu
    abd = abd / fcscal / fcscal_bd
    auu = auu / fcscal / fcscal_uu
    add = add / fcscal / fcscal_dd
    aruu = aruu / fcscal / fcscal_ruu
    ardd = ardd / fcscal / fcscal_rdd
    aluu = aluu / fcscal / fcscal_luu
    aldd = aldd / fcscal / fcscal_ldd
    afuu = afuu / fcscal / fcscal_fuu
    afdd = afdd / fcscal / fcscal_fdd
    abuu = abuu / fcscal / fcscal_buu
    abdd = abdd / fcscal / fcscal_bdd

    # Determine indices for the operator.
    ia = i - i0 + 1
    ja = j - j0 + 1
    ka = k - k0 + 1

    # Set matrix elements for bicgstab.
    ac_b[ia, ja, ka] = ac
    ar_b[ia, ja, ka] = ar
    al_b[ia, ja, ka] = al
    af_b[ia, ja, ka] = af
    ab_b[ia, ja, ka] = ab
    au_b[ia, ja, ka] = au
    ad_b[ia, ja, ka] = ad
    aru_b[ia, ja, ka] = aru
    ard_b[ia, ja, ka] = ard
    alu_b[ia, ja, ka] = alu
    ald_b[ia, ja, ka] = ald
    afu_b[ia, ja, ka] = afu
    afd_b[ia, ja, ka] = afd
    abu_b[ia, ja, ka] = abu
    abd_b[ia, ja, ka] = abd
    auu_b[ia, ja, ka] = auu
    add_b[ia, ja, ka] = add
    aruu_b[ia, ja, ka] = aruu
    ardd_b[ia, ja, ka] = ardd
    aluu_b[ia, ja, ka] = aluu
    aldd_b[ia, ja, ka] = aldd
    afuu_b[ia, ja, ka] = afuu
    afdd_b[ia, ja, ka] = afdd
    abuu_b[ia, ja, ka] = abuu
    abdd_b[ia, ja, ka] = abdd

    # Store horizontal and vertical components of AC (for
    # preconditioner).
    if preconditioner
      ach_b[ia, ja, ka] = -ar - al - af - ab
      acv_b[ia, ja, ka] = -au - ad
    end
  end

  return
end
