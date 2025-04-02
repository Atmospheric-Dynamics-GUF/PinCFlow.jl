function saturation!(state, steadystate == true)

  # do nothing 
end

function saturation!(state, dt::AbstractFloat, steadystate == false)
  # 3d extension of the saturation scheme in Boeloeni et al. (2016)

  # TODO first and fourth loop are basically identical
  (; ray, diffusion) = state.wkb
  (; sizex, sizey, sizez) = state.namelists.domain
  (; alpha_sat) = state.namelists.wkb
  (; io, jo, i0, i1, j0, j1, k0, k1, nxx, nyy, nzz) = state.domain
  (; lx, ly, lz, dx, dy, dz, jac) = state.grid
  (; rhostrattfc) = state.atmosphere

  mb2 = zeros(nxx, nyy, nzz)
  mb2k2 = zeros(nxx, nyy, nzz)

  # calculate the integrals in the numerator and denominator of 
  # K(z) in eq. (32)
  for kzrv in k0:k1
    for jyrv in j0:j1
      for ixrv in i0:i1
        if nray[ixrv, jyrv, kzrv] < 1
          continue
        end

        for iray in 1:nray[ixrv, jyrv, kzrv]

          # skip counting ray volumes with zero wave-action density
          if ray.dens[iray, ixrv, jyrv, kzrv] == 0.0
            continue
          end

          xr = ray.x[iray, ixrv, jyrv, kzrv]
          yr = ray.y[iray, ixrv, jyrv, kzrv]
          zr = ray.z[iray, ixrv, jyrv, kzrv]

          dxr = ray.dxray[iray, ixrv, jyrv, kzrv]
          dyr = ray.dyray[iray, ixrv, jyrv, kzrv]
          dzr = ray.dzray[iray, ixrv, jyrv, kzrv]

          if sizex > 1
            ix = floor((xr - lx[1]) / dx) + 1 - io
          else
            ix = 1
          end

          if sizey > 1
            jy = floor((yr - ly[1]) / dy) + 1 - jo
          else
            jy = 1
          end

          kz = kztildetfc(ix, jy, zr)

          nn_nd = stratification(state, zr, 1) # compute stratification

          wnrk = ray.k[iray, ixrv, jyrv, kzrv]
          wnrl = ray.l[iray, ixrv, jyrv, kzrv]
          wnrm = ray.m[iray, ixrv, jyrv, kzrv]

          wnrhs = wnrk^2.0 + wnrl^2.0

          dwnrk = ray.dkray[iray, ixrv, jyrv, kzrv]
          dwnrl = ray.dlray[iray, ixrv, jyrv, kzrv]
          dwnrm = ray.dmray[iray, ixrv, jyrv, kzrv]

          omir = ray.omega[iray, ixrv, jyrv, kzrv]

          densr = ray.dens[iray, ixrv, jyrv, kzrv]

          # spatial extension of ray to be taken into account
          dzi = min(dzr, jac[ix, jy, kz] * dz)
          facpsp = dzi / jac[ix, jy, kz] / dz * dwnrm

          if sizex > 1
            dxi = min(dxr, dx)
            facpsp = facpsp * dxi / dx * dwnrk
          end

          if sizey > 1
            dyi = min(dyr, dy)
            facpsp = facpsp * dyi / dy * dwnrl
          end

          integral1 = wnrhs * wnrm^2.0 / ((wnrhs + wnrm^2.0) * omir) * facpsp

          # numerator integral
          mb2[ix, jy, kz] =
            mb2[ix, jy, kz] +
            2.0 * nn_nd^2.0 / rhostrattfc[ix, jy, kz] * densr * integral1

          integral2 = wnrhs * wnrm^2.0 / omir * facpsp

          # denominator integral
          mb2k2[ix, jy, kz] =
            mb2k2[ix, jy, kz] +
            2.0 * nn_nd^2.0 / rhostrattfc[ix, jy, kz] * densr * integral2
        end
      end
    end
  end

  # calculate the turbulent eddy diffusivity K(z) in eq. (32)
  for kz in k0:k1
    for jy in j0:j1
      for ix in i0:i1
        nn_nd = stratification(state, ztfc[ix, jy, kz], 1)

        # check according to eq. (48)
        if mb2k2[ix, jy, kz] == 0.0 ||
           mb2[ix, jy, kz] < alpha_sat^2.0 * nn_nd^2.0
          diffusion[ix, jy, kz] = 0.0
        else
          diffusion[ix, jy, kz] =
            (mb2[ix, jy, kz] - alpha_sat^2.0 * nn_nd^2.0) /
            (2.0 * dt * mb2k2[ix, jy, kz])
        end
      end
    end
  end

  # loop for reducing wave action density 
  # if m^2*B^2 exceed saturation threshold
  # reduce wave action density according to eq. (51)
  for kzrv in k0:k1
    for jyrv in j0:j1
      for ixrv in j0:j1
        if nray[ixrv, jyrv, kzrv] < 1
          continue
        end

        for iray in 1:nray[ixrv, jyrv, kzrv]
          if ray.dens[iray, ixrv, jyrv, kzrv] == 0.0
            continue
          end

          xr = ray.x[iray, ixrv, jyrv, kzrv]
          yr = ray.y[iray, ixrv, jyrv, kzrv]
          zr = ray.z[iray, ixrv, jyrv, kzrv]

          if sizex > 1
            ix = floor((xr - lx[1]) / dx) + 1 - io
          else
            ix = 1
          end

          if sizey > 1
            jy = floor((yr - ly[1]) / dy) + 1 - jo
          else
            jy = 1
          end

          kz = kztildetfc(ix, jy, zr)

          wnrk = ray.k[iray, ixrv, jyrv, kzrv]
          wnrl = ray.l[iray, ixrv, jyrv, kzrv]
          wnrm = ray.m[iray, ixrv, jyrv, kzrv]

          kappa = diffusion[ix, jy, kz]

          # it can happen that the damping coefficient becomes negative
          # set it to zero in that case
          ray.dens[iray, ixrv, jyrv, kzrv] =
            ray.dens[iray, ixrv, jyrv, kzrv] *
            max(0.0, 1.0 - dt * 2.0 * kappa * 
            (wnrk^2.0 + wnrl^2.0 + wnrm^2.0))
        end
      end
    end
  end

  # diagnostic to check the impact of the saturation parameterization 
  # on the energy 
  # check m^2*B^2 again and print diagnostic 
  # exact repeat of the first loop
  mb2 = 0.0
  for kzrv in k0:k1
    for jyrv in j0:j1
      for ixrv in i0:i1
        if nray[ixrv, jyrv, kzrv] < 1
          continue
        end

        for iray in 1:nray[ixrv, jyrv, kzrv]

          # skip counting ray volumes with zero wave-action density
          if ray.dens[iray, ixrv, jyrv, kzrv] == 0.0
            continue
          end

          xr = ray.x[iray, ixrv, jyrv, kzrv]
          yr = ray.y[iray, ixrv, jyrv, kzrv]
          zr = ray.z[iray, ixrv, jyrv, kzrv]

          dxr = ray.dxray[iray, ixrv, jyrv, kzrv]
          dyr = ray.dyray[iray, ixrv, jyrv, kzrv]
          dzr = ray.dzray[iray, ixrv, jyrv, kzrv]

          if sizex > 1
            ix = floor((xr - lx[1]) / dx) + 1 - io
          else
            ix = 1
          end

          if sizey > 1
            jy = floor((yr - ly[1]) / dy) + 1 - jo
          else
            jy = 1
          end

          kz = kztildetfc(ix, jy, zr)

          nn_nd = stratification(state, zr, 1) # compute stratification

          wnrk = ray.k[iray, ixrv, jyrv, kzrv]
          wnrl = ray.l[iray, ixrv, jyrv, kzrv]
          wnrm = ray.m[iray, ixrv, jyrv, kzrv]

          wnrhs = wnrk^2.0 + wnrl^2.0

          dwnrk = ray.dkray[iray, ixrv, jyrv, kzrv]
          dwnrl = ray.dlray[iray, ixrv, jyrv, kzrv]
          dwnrm = ray.dmray[iray, ixrv, jyrv, kzrv]

          omir = ray.omega[iray, ixrv, jyrv, kzrv]

          densr = ray.dens[iray, ixrv, jyrv, kzrv]

          # spatial extension of ray to be taken into account
          dzi = min(dzr, jac[ix, jy, kz] * dz)
          facpsp = dzi / jac[ix, jy, kz] / dz * dwnrm

          if sizex > 1
            dxi = min(dxr, dx)
            facpsp = facpsp * dxi / dx * dwnrk
          end

          if sizey > 1
            dyi = min(dyr, dy)
            facpsp = facpsp * dyi / dy * dwnrl
          end

          integral1 = wnrhs * wnrm^2.0 / ((wnrhs + wnrm^2.0) * omir) * 
            facpsp

          mb2[ix, jy, kz] =
            mb2[ix, jy, kz] +
            2.0 * nn_nd^2.0 / rhostrattfc[ix, jy, kz] * densr * integral1
        end
      end
    end
  end

  for kz in k0:k1
    for jy in j0:j1
      for ix in i0:i1
        nn_nd = stratification(state, ztfc[ix, jy, kz], 1)

        if mb2[ix, jy, kz] - alpha_sat^2.0 * nn_nd^2.0 >
           1.E-3 * alpha_sat^2.0 * nn_nd^2.0
          println("SATURATION VIOLATED AT ix, jy, kz = ", ix, jy, kz)
          println("mb2(ix, jy, kz) = ", mb2[ix, jy, kz])
          println("alpha_sat ^ 2. * nn_nd ^ 2. = ", alpha_sat^2.0 * nn_nd^2.0)
        end
      end
    end
  end

  # remove ray volumes with zero wave action 
  for kz in (k0 - 1):(k1 + 1)
    for jy in (j0 - 1):(j1 + 1)
      for ix in (ix - 1):(i1 + 1)
        if nray[ix, jy, kz] < 1
          continue
        end

        nrlc = 0

        for iray in 1:nray[ix, jy, kz]
          if ray.dens[iray, ix, jy, kz] == 0.0
            continue
          end
          nrlc = nrlc + 1
          ray[nrlc, ix, jy, kz] = ray[iray, ix, jy, kz]
        end
        nray[ix, jy, kz] = nrlc
      end
    end
  end
end
