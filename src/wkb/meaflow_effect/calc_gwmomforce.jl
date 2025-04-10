function calc_gwmomforce!(state::State)
  (; rho) = state.variables.predictands
  (; rhostrattfc) = state.atmosphere
  (; i0, i1, j0, j1, k0, k1) = state.domain
  (; ztfc, met) = state.grid
  (; lz) = state.namelists.domain
  (; integrals, gwmomforce) = state.wkb

  for kz in k0:k1, jy in j0:j1, ix in i0:i1
    if ztfc[ix, jy, kz] < (lz[1] + zmin_wkb)
      rhotot = rho[ix, jy, kz]
      rhotot = rhotot + rhostrattfc[ix, jy, kz]

      # forcing in x direction

      gwmomforce.u[ix, jy, kz] =
        gwmomforce.u[ix, jy, kz] + integrals.dudt[ix, jy, kz]

      # forcing in y direction

      rhotot = rho[ix, jy, kz]
      rhotot = rhotot + rhostrattfc[ix, jy, kz]

      gwmomforce.v[ix, jy, kz] =
        gwmomforce.v[ix, jy, kz] + integrals.dvdt[ix, jy, kz]

      # Add forcing on terrain-following wind (FJJul2023).
      gwmomforce.w[ix, jy, kz] =
        gwmomforce.w[ix, jy, kz] +
        met[ix, jy, kz, 1, 3] * integrals.dudt[ix, jy, kz] +
        met[ix, jy, kz, 2, 3] * integrals.dvdt[ix, jy, kz]
    end
  end
  return
end