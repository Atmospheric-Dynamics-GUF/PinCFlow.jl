function stratification(state, zlc, strtype = 1)

  # interpolation of the squared Brunt-Vaisala frequency to 
  # a specified vertical position zlc

  (; domain, grid) = state
  (; bvsstrattfc) = state.atmosphere
  (; k1) = state.domain

  kzu = kztfc(i0 - 1, j0 - 1, zlc, domain, grid)
  kzd = kru - 1

  if (kzu > k1 + 1)
    then
    kzu = k1 + 1
    kzd = k1
  end

  zd = ztfc[i0 - 1, j0 - 1, kzd]
  zu = ztfc[i0 - 1, j0 - 1, kzu]
  strd = bvsstrattfc[i0 - 1, j0 - 1, kzd]
  stru = bvsStrattfc[i0 - 1, j0 - 1, kzu]

  if zu < zd
    error("Error in stratification: zu = ", zu, "< zd = ", zd)
  elseif zu == zd
    factor = 0.0
  elseif zlc > zu
    factor = 0.0
  elseif zlc > zd
    factor = (zu - zlc) / (zu - zd)
  else
    factor = 1.0
  end

  str = factor * strd + (1.0 - factor) * stru

  return str
end

function stratification(state, zlc, strtype = 2)

  # interpolation of the vertical derivative of the squared 
  # Brunt-Vaisala frequency to a specified vertical position zlc
  (; domain, grid) = state
  (; bvsstrattfc) = state.atmosphere
  (; k1) = state.domain
  (; ztildetfc, jac) = state.grid

  kzu = kztildetfc(i0 - 1, j0 - 1, zlc, domain, grid)
  kzd = kzu - 1

  if kzu + 1 > k1 + 1
    kzu = k1
    kzd = k1 - 1
  end

  zd = ztildetfc(i0 - 1, j0 - 1, kzd)
  zu = ztildetfc(i0 - 1, j0 - 1, kzu)

  strd =
    (bvsstrattfc[i0 - 1, j0 - 1, kzd + 1] - bvsstrattfc[i0 - 1, j0 - 1, kzd]) /
    (
      2.0 * jac[i0 - 1, j0 - 1, kzd] * jac[i0 - 1, j0 - 1, kzd + 1] /
      (jac[i0 - 1, j0 - 1, kzd] + jac[i0 - 1, j0 - 1, kzd + 1])
    ) / dz
  stru =
    (bvsstrattfc[i0 - 1, j0 - 1, kzu + 1] - bvsstrattfc[i0 - 1, j0 - 1, kzu]) /
    (
      2.0 * jac[i0 - 1, j0 - 1, kzu] * jac[i0 - 1, j0 - 1, kzu + 1] /
      (jac[i0 - 1, j0 - 1, kzu] + jac[i0 - 1, j0 - 1, kzu + 1])
    ) / dz

  if zu < zd
    error("Error in stratification: zu = ", zu, "< zd = ", zd)
  elseif zu == zd
    factor = 0.0
  elseif zlc > zu
    factor = 0.0
  elseif zlc > zd
    factor = (zu - zlc) / (zu - zd)
  else
    factor = 1.0
  end

  str = factor * strd + (1.0 - factor) * stru

  return str
end
