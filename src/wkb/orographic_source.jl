function orographic_source!(
  namelists::Namelists,
  constants::Constants,
  domain::Domain,
  grid::Grid,
  atmosphere::Atmosphere,
  predictands::Predictands,
  omi_ini::AbstractArray{<:AbstractFloat, 4},
  wnk_ini::AbstractArray{<:AbstractFloat, 4},
  wnl_ini::AbstractArray{<:AbstractFloat, 4},
  wnm_ini::AbstractArray{<:AbstractFloat, 4},
  wad_ini::AbstractArray{<:AbstractFloat, 4},
  zb::AbstractMatrix{<:AbstractFloat},
)

  # Get all necessary fields.
  (; f_cor_nd) = namelists.atmosphere
  (; branchr, blocking, long_threshold) = namelists.wkb
  (; tref) = constants
  (; i0, i1, j0, j1, k0, k1) = domain
  (; dz, jac, ztildetfc, k_spectrum, l_spectrum, topography_spectrum) = grid
  (; rhostrattfc, bvsstrattfc) = atmosphere
  (; u, v) = predictands

  # Set Coriolis parameter.
  f_cor_nd = f_coriolis_dim * tref

  # Iterate over surface grid cells.
  for jy in j0:j1, ix in i0:i1

    # Average mean wind, reference density and buoyancy frequency. This should
    # be done without a vertical loop.
    uavg = 0.0
    vavg = 0.0
    rhoavg = 0.0
    bvsavg = 0.0
    dzsum = 0.0
    for kz in k0:k1
      uavg += 0.5 * (u[ix, jy, kz] + u[ix - 1, jy, kz]) * jac[ix, jy, kz] * dz
      vavg += 0.5 * (v[ix, jy, kz] + v[ix, jy - 1, kz]) * jac[ix, jy, kz] * dz
      rhoavg += rhostrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
      bvsavg += bvsstrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
      dzsum += jac[ix, jy, kz] * dz
      if ztildetfc[ix, jy, kz] >
         ztildetfc[ix, jy, k0 - 1] + sum(abs.(topography_spectrum[:, ix, jy]))
        break
      end
    end
    uavg = uavg / dzsum
    vavg = vavg / dzsum
    rhoavg = rhoavg / dzsum
    bvsavg = bvsavg / dzsum

    # Determine the blocked layer.
    if blocking && sum(abs(topography_spectrum[:, ix, jy])) > 0.0
      long =
        sqrt(bvsavg) / sqrt(uavg^2.0 + vavg^2.0) *
        sum(abs.(topography_spectrum[:, ix, jy]))
      ratio = min(1.0, long_threshold / long)
      zb[ix, jy] =
        ztildetfc[ix, jy, 0] +
        sum(abs.(topography_spectrum[:, ix, jy])) * (1.0 - 2.0 * ratio)
    elseif blocking
      ratio = 1.0
      zb[ix, jy] = ztildetfc[ix, jy, 0]
    else
      ratio = 1.0
    end

    # Set launch level.
    kz = k0 - 1

    # Iterate over wave modes.
    for iwm in 1:nwm

      # Set wavenumbers.
      wnk = k_spectrum[iwm, ix, jy]
      wnl = l_spectrum[iwm, ix, jy]
      wnh = sqrt(wnk^2.0 + wnl^2.0)

      # Compute intrinsic frequency from orographic wavenumbers.
      omi = -uavg * wnk - vavg * wnl

      # Adjust the signs to be consistent with the chosen frequency
      # branch.
      if omi * branchr < 0.0
        omi = -omi
        wnk = -wnk
        wnl = -wnl
      end

      # Compute vertical wavenumber and wave-action density.
      if omi^2 > f_cor_nd^2 && omi^2 < bvsavg

        # Compute vertical wavenumber.
        wnm = -branchr * sqrt(wnh^2 * (bvsavg - omi^2) / (omi^2 - f_cor_nd^2))

        # Compute displacement.
        displm = ratio * abs(topography_spectrum[iwm, ix, jy])

        # Compute wave-action density.
        wad = 0.5 * rhoavg * displm^2 * omi * (wnh^2 + wnm^2) / wnh^2

        # Set to zero if something went wrong.
        if wad != wad || wnm != wnm
          wad = 0.0
          wnm = 0.0
        end

        # Account for critical and reflecting levels.
      else
        wad = 0.0
        wnm = 0.0
      end

      # Save intrinsic frequency, wavenumbers and wave-action density.
      omi_ini[iwm, ix, jy, kz] = omi
      wnk_ini[iwm, ix, jy, kz] = wnk
      wnl_ini[iwm, ix, jy, kz] = wnl
      wnm_ini[iwm, ix, jy, kz] = wnm
      wad_ini[iwm, ix, jy, kz] = wad
    end
  end
end
