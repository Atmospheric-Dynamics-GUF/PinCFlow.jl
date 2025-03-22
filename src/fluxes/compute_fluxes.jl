function compute_fluxes!(state::State, predictands::Predictands)
  compute_fluxes!(state, predictands, Rho())
  compute_fluxes!(state, predictands, RhoP())
  compute_fluxes!(state, predictands, U())
  compute_fluxes!(state, predictands, V())
  compute_fluxes!(state, predictands, W())
  return
end

function compute_fluxes!(state::State, predictands::Predictands, variable::Rho)

  # Get all necessary fields.
  (; domain, atmosphere) = state
  (; nx, ny, nz) = domain
  (; jac) = state.grid
  (; pstrattfc, rhostrattfc) = atmosphere
  (; rhotilde) = state.variables.reconstructions
  (; phirho) = state.variables.fluxes

  # Get old wind.
  (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

  #-----------------------------------------
  #             Zonal fluxes
  #-----------------------------------------

  for k in 1:nz, j in 1:ny, i in 0:nx
    rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
    pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
    rhor = rhotilde[i + 1, j, k, 1, 0] + rhostratedger / pedger
    rhol = rhotilde[i, j, k, 1, 1] + rhostratedger / pedger

    pedger =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
      )
    usurf = pedger * u0[i, j, k]

    frho = compute_flux(usurf, rhol, rhor)

    phirho[i, j, k, 1] = frho
  end

  #-----------------------------------------
  #           Meridional fluxes
  #-----------------------------------------

  for k in 1:nz, j in 0:ny, i in 1:nx
    rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
    pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
    rhof = rhotilde[i, j + 1, k, 2, 0] + rhostratedgef / pedgef
    rhob = rhotilde[i, j, k, 2, 1] + rhostratedgef / pedgef

    pedgef =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
      )
    vsurf = pedgef * v0[i, j, k]

    grho = compute_flux(vsurf, rhob, rhof)

    phirho[i, j, k, 2] = grho
  end

  #-----------------------------------------
  #            Vertical fluxes
  #-----------------------------------------

  for k in 0:nz, j in 1:ny, i in 1:nx
    rhostratedgeu =
      (
        jac[i, j, k + 1] * rhostrattfc[i, j, k] +
        jac[i, j, k] * rhostrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    pedgeu =
      (
        jac[i, j, k + 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    rhou = rhotilde[i, j, k + 1, 3, 0] + rhostratedgeu / pedgeu
    rhod = rhotilde[i, j, k, 3, 1] + rhostratedgeu / pedgeu

    pedgeu =
      jac[i, j, k] *
      jac[i, j, k + 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1])
    wsurf = pedgeu * w0[i, j, k]

    hrho = compute_flux(wsurf, rhod, rhou)

    phirho[i, j, k, 3] = hrho
  end

  # Return.
  return
end

function compute_fluxes!(state::State, predictands::Predictands, variable::RhoP)

  # Get all necessary fields.
  (; domain, atmosphere) = state
  (; nx, ny, nz) = domain
  (; jac) = state.grid
  (; pstrattfc, rhostrattfc) = atmosphere
  (; rhoptilde) = state.variables.reconstructions
  (; phirhop) = state.variables.fluxes

  # Get old wind.
  (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

  #-----------------------------------------
  #             Zonal fluxes
  #-----------------------------------------

  for k in 1:nz, j in 1:ny, i in 0:nx
    rhor = rhoptilde[i + 1, j, k, 1, 0]
    rhol = rhoptilde[i, j, k, 1, 1]

    pedger =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
      )
    usurf = pedger * u0[i, j, k]

    frhop = compute_flux(usurf, rhol, rhor)

    phirhop[i, j, k, 1] = frhop
  end

  #-----------------------------------------
  #           Meridional fluxes
  #-----------------------------------------

  for k in 1:nz, j in 0:ny, i in 1:nx
    rhof = rhoptilde[i, j + 1, k, 2, 0]
    rhob = rhoptilde[i, j, k, 2, 1]

    pedgef =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
      )
    vsurf = pedgef * v0[i, j, k]

    grhop = compute_flux(vsurf, rhob, rhof)

    phirhop[i, j, k, 2] = grhop
  end

  #-----------------------------------------
  #            Vertical fluxes
  #-----------------------------------------

  for k in 0:nz, j in 1:ny, i in 1:nx
    rhou = rhoptilde[i, j, k + 1, 3, 0]
    rhod = rhoptilde[i, j, k, 3, 1]

    pedgeu =
      jac[i, j, k] *
      jac[i, j, k + 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1])
    wsurf = pedgeu * w0[i, j, k]

    hrhop = compute_flux(wsurf, rhod, rhou)

    phirhop[i, j, k, 3] = hrhop
  end

  # Return.
  return
end

function compute_fluxes!(
  state::State,
  old_predictands::Predictands,
  variable::U,
)

  # Get all necessary fields.
  (; re) = state.constants
  (; nx, ny, nz) = state.domain
  (; dx, dy, dz, jac, met) = state.grid
  (; pstrattfc, rhostrattfc) = state.atmosphere
  (; utilde) = state.variables.reconstructions
  (; phiu) = state.variables.fluxes
  (; predictands) = state.variables

  # Get old wind.
  (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

  #-----------------------------------------
  #             Zonal fluxes
  #-----------------------------------------

  for k in 1:nz, j in 1:ny, i in (-1):nx
    # The uTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    ur = utilde[i + 1, j, k, 1, 0]
    ul = utilde[i, j, k, 1, 1]

    pedger =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
      )
    predger =
      0.5 * (
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k] +
        jac[i + 2, j, k] * pstrattfc[i + 2, j, k]
      )
    usurf = 0.5 * (pedger * u0[i, j, k] + predger * u0[i + 1, j, k])

    frhou = compute_flux(usurf, ul, ur)

    phiu[i, j, k, 1] = frhou
  end

  #-----------------------------------------
  #           Meridional fluxes
  #-----------------------------------------

  for k in 1:nz, j in 0:ny, i in 0:nx
    # The uTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    uf = utilde[i, j + 1, k, 2, 0]
    ub = utilde[i, j, k, 2, 1]

    pedgef =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
      )
    predgef =
      0.5 * (
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k] +
        jac[i + 1, j + 1, k] * pstrattfc[i + 1, j + 1, k]
      )
    vsurf = 0.5 * (pedgef * v0[i, j, k] + predgef * v0[i + 1, j, k])

    grhou = compute_flux(vsurf, ub, uf)

    phiu[i, j, k, 2] = grhou
  end

  #-----------------------------------------
  #            Vertical fluxes
  #-----------------------------------------

  for k in 0:nz, j in 1:ny, i in 0:nx
    # The uTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    uu = utilde[i, j, k + 1, 3, 0]
    ud = utilde[i, j, k, 3, 1]

    pedgeu =
      jac[i, j, k] *
      jac[i, j, k + 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1])
    predgeu =
      jac[i + 1, j, k] *
      jac[i + 1, j, k + 1] *
      (pstrattfc[i + 1, j, k] + pstrattfc[i + 1, j, k + 1]) /
      (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
    wsurf = 0.5 * (pedgeu * w0[i, j, k] + predgeu * w0[i + 1, j, k])

    hrhou = compute_flux(wsurf, ud, uu)

    phiu[i, j, k, 3] = hrhou
  end

  #-------------------------------------------------------------------
  #                          Viscous fluxes
  #-------------------------------------------------------------------

  if re >= 1.0e20
    return
  end

  #-----------------------------------------
  #             Zonal fluxes
  #-----------------------------------------

  for k in 1:nz, j in 1:ny, i in (-1):nx
    coef_v = 1 / re * rhostrattfc[i + 1, j, 1]

    frhou_visc =
      coef_v *
      jac[i + 1, j, k] *
      compute_stress_tensor(i + 1, j, k, 1, 1, predictands, grid)

    phiu[i, j, k, 1] -= frhou_visc
  end

  #-----------------------------------------
  #           Meridional fluxes
  #-----------------------------------------

  for k in 1:nz, j in 0:ny, i in 0:nx
    coef_v =
      1 / re *
      0.25 *
      (
        rhostrattfc[i, j, 1] +
        rhostrattfc[i + 1, j, 1] +
        rhostrattfc[i, j + 1, 1] +
        rhostrattfc[i + 1, j + 1, 1]
      )

    grhou_visc =
      coef_v *
      0.25 *
      (
        jac[i, j, k] * compute_stress_tensor(i, j, k, 1, 2, predictands, grid) +
        jac[i + 1, j, k] *
        compute_stress_tensor(i + 1, j, k, 1, 2, predictands, grid) +
        jac[i, j + 1, k] *
        compute_stress_tensor(i, j + 1, k, 1, 2, predictands, grid) +
        jac[i + 1, j + 1, k] *
        compute_stress_tensor(i + 1, j + 1, k, 1, 2, predictands, grid)
      )

    phiu[i, j, k, 2] -= grhou_visc
  end

  #-----------------------------------------
  #            Vertical fluxes
  #-----------------------------------------

  for k in 0:nz, j in 1:ny, i in 0:nx
    coef_v = 1 / re * 0.5 * (rhostrattfc[i, j, 1] + rhostrattfc[i + 1, j, 1])

    stresstens13 =
      met[i, j, k, 1, 3] *
      compute_stress_tensor(i, j, k, 1, 1, predictands, grid) +
      met[i, j, k, 2, 3] *
      compute_stress_tensor(i, j, k, 1, 2, predictands, grid) +
      compute_stress_tensor(i, j, k, 1, 3, predictands, grid) / jac[i, j, k]
    stresstens13r =
      met[i + 1, j, k, 1, 3] *
      compute_stress_tensor(i + 1, j, k, 1, 1, predictands, grid) +
      met[i + 1, j, k, 2, 3] *
      compute_stress_tensor(i + 1, j, k, 1, 2, predictands, grid) +
      compute_stress_tensor(i + 1, j, k, 1, 3, predictands, grid) /
      jac[i + 1, j, k]
    stresstens13u =
      met[i, j, k + 1, 1, 3] *
      compute_stress_tensor(i, j, k + 1, 1, 1, predictands, grid) +
      met[i, j, k + 1, 2, 3] *
      compute_stress_tensor(i, j, k + 1, 1, 2, predictands, grid) +
      compute_stress_tensor(i, j, k + 1, 1, 3, predictands, grid) /
      jac[i, j, k + 1]
    stresstens13ru =
      met[i + 1, j, k + 1, 1, 3] *
      compute_stress_tensor(i + 1, j, k + 1, 1, 1, predictands, grid) +
      met[i + 1, j, k + 1, 2, 3] *
      compute_stress_tensor(i + 1, j, k + 1, 1, 2, predictands, grid) +
      compute_stress_tensor(i + 1, j, k + 1, 1, 3, predictands, grid) /
      jac[i + 1, j, k + 1]
    hrhou_visc =
      coef_v *
      0.5 *
      (
        jac[i, j, k] * jac[i, j, k + 1] * (stresstens13 + stresstens13u) /
        (jac[i, j, k] + jac[i, j, k + 1]) +
        jac[i + 1, j, k] *
        jac[i + 1, j, k + 1] *
        (stresstens13r + stresstens13ru) /
        (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
      )

    phiu[i, j, k, 3] -= hrhou_visc
  end

  # Return.
  return
end

function compute_fluxes!(
  state::State,
  old_predictands::Predictands,
  variable::V,
)

  # Get all necessary fields.
  (; re) = state.constants
  (; nx, ny, nz) = state.domain
  (; dx, dy, dz, jac, met) = state.grid
  (; pstrattfc, rhostrattfc) = state.atmosphere
  (; vtilde) = state.variables.reconstructions
  (; phiv) = state.variables.fluxes
  (; predictands) = state.variables

  # Get old wind.
  (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

  #-----------------------------------------
  #             Zonal fluxes
  #-----------------------------------------

  for k in 1:nz, j in 0:ny, i in 0:nx
    # The vTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    vr = vtilde[i + 1, j, k, 1, 0]
    vl = vtilde[i, j, k, 1, 1]

    pedger =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
      )
    pfedger =
      0.5 * (
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k] +
        jac[i + 1, j + 1, k] * pstrattfc[i + 1, j + 1, k]
      )
    usurf = 0.5 * (pedger * u0[i, j, k] + pfedger * u0[i, j + 1, k])

    frhov = compute_flux(usurf, vl, vr)

    phiv[i, j, k, 1] = frhov
  end

  #-----------------------------------------
  #           Meridional fluxes
  #-----------------------------------------

  for k in 1:nz, j in (-1):ny, i in 1:nx
    # The vTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    vf = vtilde[i, j + 1, k, 2, 0]
    vb = vtilde[i, j, k, 2, 1]

    pedgef =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
      )
    pfedgef =
      0.5 * (
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k] +
        jac[i, j + 2, k] * pstrattfc[i, j + 2, k]
      )
    vsurf = 0.5 * (pedgef * v0[i, j, k] + pfedgef * v0[i, j + 1, k])

    grhov = compute_flux(vsurf, vb, vf)

    phiv[i, j, k, 2] = grhov
  end

  #-----------------------------------------
  #            Vertical fluxes
  #-----------------------------------------

  for k in 0:nz, j in 0:ny, i in 1:nx
    # The vTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    vu = vtilde[i, j, k + 1, 3, 0]
    vd = vtilde[i, j, k, 3, 1]

    pedgeu =
      jac[i, j, k] *
      jac[i, j, k + 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1])
    pfedgeu =
      jac[i, j + 1, k] *
      jac[i, j + 1, k + 1] *
      (pstrattfc[i, j + 1, k] + pstrattfc[i, j + 1, k + 1]) /
      (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
    wsurf = 0.5 * (pedgeu * w0[i, j, k] + pfedgeu * w0[i, j + 1, k])

    hrhov = compute_flux(wsurf, vd, vu)

    phiv[i, j, k, 3] = hrhov
  end

  #-------------------------------------------------------------------
  #                          Viscous fluxes
  #-------------------------------------------------------------------

  if re >= 1.0e20
    return
  end

  #-----------------------------------------
  #             Zonal fluxes
  #-----------------------------------------

  for k in 1:nz, j in 0:ny, i in 0:nx
    coef_v =
      1 / re *
      0.25 *
      (
        rhostrattfc[i, j, 1] +
        rhostrattfc[i + 1, j, 1] +
        rhostrattfc[i, j + 1, 1] +
        rhostrattfc[i + 1, j + 1, 1]
      )

    frhov_visc =
      coef_v *
      0.25 *
      (
        jac[i, j, k] * compute_stress_tensor(i, j, k, 2, 1, predictands, grid) +
        jac[i + 1, j, k] *
        compute_stress_tensor(i + 1, j, k, 2, 1, predictands, grid) +
        jac[i, j + 1, k] *
        compute_stress_tensor(i, j + 1, k, 2, 1, predictands, grid) +
        jac[i + 1, j + 1, k] *
        compute_stress_tensor(i + 1, j + 1, k, 2, 1, predictands, grid)
      )

    phiv[i, j, k, 1] -= frhov_visc
  end

  #-----------------------------------------
  #           Meridional fluxes
  #-----------------------------------------

  for k in 1:nz, j in (-1):ny, i in 1:nx
    coef_v = 1 / re * rhostrattfc[i, j + 1, 1]

    grhov_visc =
      coef_v *
      jac[i, j + 1, k] *
      compute_stress_tensor(i, j + 1, k, 2, 2, predictands, grid)

    phiv[i, j, k, 2] -= grhov_visc
  end

  #-----------------------------------------
  #            Vertical fluxes
  #-----------------------------------------

  for k in 0:nz, j in 0:ny, i in 1:nx
    coef_v = 1 / re * 0.5 * (rhostrattfc[i, j, 1] + rhostrattfc[i, j + 1, 1])

    stresstens23 =
      met[i, j, k, 1, 3] *
      compute_stress_tensor(i, j, k, 2, 1, predictands, grid) +
      met[i, j, k, 2, 3] *
      compute_stress_tensor(i, j, k, 2, 2, predictands, grid) +
      compute_stress_tensor(i, j, k, 2, 3, predictands, grid) / jac[i, j, k]
    stresstens23f =
      met[i, j + 1, k, 1, 3] *
      compute_stress_tensor(i, j + 1, k, 2, 1, predictands, grid) +
      met[i, j + 1, k, 2, 3] *
      compute_stress_tensor(i, j + 1, k, 2, 2, predictands, grid) +
      compute_stress_tensor(i, j + 1, k, 2, 3, predictands, grid) /
      jac[i, j + 1, k]
    stresstens23u =
      met[i, j, k + 1, 1, 3] *
      compute_stress_tensor(i, j, k + 1, 2, 1, predictands, grid) +
      met[i, j, k + 1, 2, 3] *
      compute_stress_tensor(i, j, k + 1, 2, 2, predictands, grid) +
      compute_stress_tensor(i, j, k + 1, 2, 3, predictands, grid) /
      jac[i, j, k + 1]
    stresstens23fu =
      met[i, j + 1, k + 1, 1, 3] *
      compute_stress_tensor(i, j + 1, k + 1, 2, 1, predictands, grid) +
      met[i, j + 1, k + 1, 2, 3] *
      compute_stress_tensor(i, j + 1, k + 1, 2, 2, predictands, grid) +
      compute_stress_tensor(i, j + 1, k + 1, 2, 3, predictands, grid) /
      jac[i, j + 1, k + 1]
    hrhov_visc =
      coef_v *
      0.5 *
      (
        jac[i, j, k] * jac[i, j, k + 1] * (stresstens23 + stresstens23u) /
        (jac[i, j, k] + jac[i, j, k + 1]) +
        jac[i, j + 1, k] *
        jac[i, j + 1, k + 1] *
        (stresstens23f + stresstens23fu) /
        (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
      )

    phiv[i, j, k, 3] -= hrhov_visc
  end

  # Return.
  return
end

function compute_fluxes!(
  state::State,
  old_predictands::Predictands,
  variable::W,
)

  # Get all necessary fields.
  (; re) = state.constants
  (; nx, ny, nz) = state.domain
  (; dx, dy, dz, jac, met) = state.grid
  (; pstrattfc, rhostrattfc) = state.atmosphere
  (; wtilde) = state.variables.reconstructions
  (; phiw) = state.variables.fluxes
  (; predictands) = state.variables

  # Get old wind.
  (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

  #-----------------------------------------
  #             Zonal fluxes
  #-----------------------------------------

  for k in 0:nz, j in 1:ny, i in 0:nx
    # The wTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    wr = wtilde[i + 1, j, k, 1, 0]
    wl = wtilde[i, j, k, 1, 1]

    pedger =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i + 1, j, k] * pstrattfc[i + 1, j, k]
      )
    puedger =
      0.5 * (
        jac[i, j, k + 1] * pstrattfc[i, j, k + 1] +
        jac[i + 1, j, k + 1] * pstrattfc[i + 1, j, k + 1]
      )
    usurf =
      (
        (jac[i, j, k + 1] + jac[i + 1, j, k + 1]) * pedger * u0[i, j, k] +
        (jac[i, j, k] + jac[i + 1, j, k]) * puedger * u0[i, j, k + 1]
      ) / (
        jac[i, j, k] +
        jac[i + 1, j, k] +
        jac[i, j, k + 1] +
        jac[i + 1, j, k + 1]
      )

    frhow = compute_flux(usurf, wl, wr)

    phiw[i, j, k, 1] = frhow
  end

  #-----------------------------------------
  #           Meridional fluxes
  #-----------------------------------------

  for k in 0:nz, j in 0:ny, i in 1:nx
    # The wTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    wf = wtilde[i, j + 1, k, 2, 0]
    wb = wtilde[i, j, k, 2, 1]

    pedgef =
      0.5 * (
        jac[i, j, k] * pstrattfc[i, j, k] +
        jac[i, j + 1, k] * pstrattfc[i, j + 1, k]
      )
    puedgef =
      0.5 * (
        jac[i, j, k + 1] * pstrattfc[i, j, k + 1] +
        jac[i, j + 1, k + 1] * pstrattfc[i, j + 1, k + 1]
      )
    vsurf =
      (
        (jac[i, j, k + 1] + jac[i, j + 1, k + 1]) * pedgef * v0[i, j, k] +
        (jac[i, j, k] + jac[i, j + 1, k]) * puedgef * v0[i, j, k + 1]
      ) / (
        jac[i, j, k] +
        jac[i, j + 1, k] +
        jac[i, j, k + 1] +
        jac[i, j + 1, k + 1]
      )

    grhow = compute_flux(vsurf, wb, wf)

    phiw[i, j, k, 2] = grhow
  end

  #-----------------------------------------
  #            Vertical fluxes
  #-----------------------------------------

  for k in (-1):nz, j in 1:ny, i in 1:nx
    # The wTilde are the reconstructed specific momenta, divided by P.
    # These are to be multiplied by the linearly interpolated velocities
    # (times P) in order to obtain the desired momentum fluxes.

    wu = wtilde[i, j, k + 1, 3, 0]
    wd = wtilde[i, j, k, 3, 1]

    pedgeu =
      jac[i, j, k] *
      jac[i, j, k + 1] *
      (pstrattfc[i, j, k] + pstrattfc[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1])
    puedgeu =
      jac[i, j, k + 1] *
      jac[i, j, k + 2] *
      (pstrattfc[i, j, k + 1] + pstrattfc[i, j, k + 2]) /
      (jac[i, j, k + 1] + jac[i, j, k + 2])
    wsurf = 0.5 * (pedgeu * w0[i, j, k] + puedgeu * w0[i, j, k + 1])

    hrhow = compute_flux(wsurf, wd, wu)

    phiw[i, j, k, 3] = hrhow
  end

  #-------------------------------------------------------------------
  #                          Viscous fluxes
  #-------------------------------------------------------------------

  if re >= 1.0e20
    return
  end

  #-----------------------------------------
  #             Zonal fluxes
  #-----------------------------------------

  for k in 0:nz, j in 1:ny, i in 0:nx
    coef_v = 1 / re * 0.5 * (rhostrattfc[i, j, 1] + rhostrattfc[i + 1, j, 1])

    frhow_visc =
      coef_v *
      0.5 *
      (
        jac[i, j, k] *
        jac[i, j, k + 1] *
        (
          compute_stress_tensor(i, j, k, 3, 1, predictands, grid) +
          compute_stress_tensor(i, j, k + 1, 3, 1, predictands, grid)
        ) / (jac[i, j, k] + jac[i, j, k + 1]) +
        jac[i + 1, j, k] *
        jac[i + 1, j, k + 1] *
        (
          compute_stress_tensor(i + 1, j, k, 3, 1, predictands, grid) +
          compute_stress_tensor(i + 1, j, k + 1, 3, 1, predictands, grid)
        ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
      )

    phiw[i, j, k, 1] -= frhow_visc
  end

  #-----------------------------------------
  #           Meridional fluxes
  #-----------------------------------------

  for k in 0:nz, j in 0:ny, i in 1:nx
    coef_v = 1 / re * 0.5 * (rhostrattfc[i, j, 1] + rhostrattfc[i, j + 1, 1])

    grhow_visc =
      coef_v *
      0.5 *
      (
        jac[i, j, k] *
        jac[i, j, k + 1] *
        (
          compute_stress_tensor(i, j, k, 3, 1, predictands, grid) +
          compute_stress_tensor(i, j, k + 1, 3, 1, predictands, grid)
        ) / (jac[i, j, k] + jac[i, j, k + 1]) +
        jac[i, j + 1, k] *
        jac[i, j + 1, k + 1] *
        (
          compute_stress_tensor(i, j + 1, k, 3, 1, predictands, grid) +
          compute_stress_tensor(i, j + 1, k + 1, 3, 1, predictands, grid)
        ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
      )

    phiw[i, j, k, 2] -= grhow_visc
  end

  #-----------------------------------------
  #            Vertical fluxes
  #-----------------------------------------

  for k in (-1):nz, j in 1:ny, i in 1:nx
    coef_v = 1 / re * rhostrattfc[i, j, 1]

    hrhow_visc =
      coef_v * (
        jac[i, j, k + 1] *
        met[i, j, k + 1, 1, 3] *
        compute_stress_tensor(i, j, k + 1, 3, 1, predictands, grid) +
        jac[i, j, k + 1] *
        met[i, j, k + 1, 2, 3] *
        compute_stress_tensor(i, j, k + 1, 3, 2, predictands, grid) +
        compute_stress_tensor(i, j, k + 1, 3, 3, predictands, grid)
      )

    phiw[i, j, k, 3] -= hrhow_visc
  end

  # Return.
  return
end
