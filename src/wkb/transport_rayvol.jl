abstract type AbstractDir end
struct xDir <: AbstractDir end
struct yDir <: AbstractDir end
struct zDir <: AbstractDir end

function as_index(dir::AbstractDir)
if dir == xDir()
    return 1
elseif dir == yDir()
    return 2
elseif dir == zDir()
    return 3
end
end

function size(domain, dir)
    if dir == xDir()
        return domain.sizex
    elseif dir == yDir()
        return domain.sizey
    elseif dir == zDir()
        return domain.sizez
    end
end



function interp_meanflow(xr, yr, zr, var, dir::xDir)
    cgrx1 = cgirx + meanflow(xr1, yr, zr, var, 1)
    cgrx2 = cgirx + meanflow(xr2, yr, zr, var, 1)
    # ! group velocity in x direction for the carrier ray
    cgrx = 0.5 * (cgrx1 + cgrx2)
    return cgrx
end

function interp_meanflow(dir::yDir)
    meanflow(xr, yr1, zr, var, 2, vyr1)
    meanflow(xr, yr2, zr, var, 2, vyr2)
    # ! group velocity in y direction at the two edges in y
    cgry1 = cgiry + vyr1
    cgry2 = cgiry + vyr2
    # ! group velocity in y direction for the carrier ray
    cgry = 0.5 * (cgry1 + cgry2)
    return cgry
end

function interp_meanflow(dir::zDir)
    cgrz1 = cgirz1 # ! +wzr1
    cgrz2 = cgirz2 # ! +wzr2
    # ! group velocity in z direction for the carrier ray
    cgrz = 0.5 * (cgrz1 + cgrz2)
    return cgrz
end

function displacement(rijk, max_groupvs, rays, domain, dir::AbstractDir)
 if domainsize(domain, dir) > 1 && rijk[4] > 0 && !single_column
     dir_idx = dir_to_int(dir)
     cgr = interp_meanflow(dir)
     F = cgr # !allow horizontal ray propagation
     dxRay[dir_idx, rijk] = dt * f + alphark[rkstage] * dxray[dir_idx, rijk]
     rays.x[rijk] += betark[rkstage] * dxray[dir_idx, rijk]
     # ! update maximum group velocity in x direction
     max_groupvs = max(max_group_vs, abs(cgr))
 end
 return max_groupvs
end

function displacement(rijk, max_group_vs, cgz_max_tfc, rays, domain, ::zDir)

    # !-----------------------------
    # !     vertical displacement
    # !-----------------------------

    # ! RK update

    # ! in line with the asymptotic results the vertcal wind is
    # ! NOT added to the intrinsic vertical group velocity
    # ! should one want to change this, one would also have to
    # ! take the vertical-wind gradient into account in the
    # ! prognostic equations for the wave number

    # ! call meanflow(xr,yr,zr1,var,3,wzr1)
    # ! call meanflow(xr,yr,zr2,var,3,wzr2)

    # ! group velocity in z direction at the two edges in z
    cgrz1 = cgirz1 # ! +wzr1
    cgrz2 = cgirz2 # ! +wzr2
    # ! group velocity in z direction for the carrier ray
    cgrz = 0.5 * (cgrz1 + cgrz2)
    F = cgrz
    dxRay[3, rijk] = dt * F + alphark[rkstage] * dxray[3, rijk]
    rays.z[rijk] += betark[rkstage] * dxray[3, rijk]

    # ! update maximum group velocity in z direction
    max_group_vs = max(max_group_vs, abs(cgrz))
    ijk = CartesianIndex(rijk[2], rijk[3], rijk[4])
    cgz_max_tfc[ijk] = max(cgz_max_tfc[ijk], abs(cgrz))
    return max_groupvs
end


function displacement(rijk, rays, max_group_vs, domain)
    cgx_max = displacement(rijk, rays, domain, max_group_vs[1], xDir())
    cgy_max= displacement(rijk, rays, domain, max_group_vs[2], yDir())
    cgz_max = displacement(rijk, rays, domain, max_group_vs[3], cgz_max_tfc, zDir())
    return (cgx_max, cgy_max, cgz_max)
end

function intrinsic_group_velocity(domain)
    sizex, sizey, sizez = size(domain)

    if (sizex > 1)
      # ! intrinsic group velocity in x direction not depending
      # ! on x
      cgirx = wnrk * (nnr - omir^2) / (omir * (wnrh^2 + wnrm^2))
    end

    if (sizey > 1)
      # ! intrinsic group velocity in y direction not depending
      # ! on y
      cgiry = wnrl * (nnr - omir^2) / (omir * (wnrh^2 + wnrm^2))
    end

    cgirz1 = -wnrm * (omir1^2 - f_cor_nd^2) / (omir1 * (wnrh^2 + wnrm^2))
    cgirz2 = -wnrm * (omir2^2 - f_cor_nd^2) / (omir2 * (wnrh^2 + wnrm^2))

    return cgrz
end



function transport_unsteady(rkstage, state)
  # ! initialize RK-tendencies at first RK stage
  (; case_wkb)=  state.namelists.wkb
  (; nx,ny,nz) = state.domain
  (; nray) = state.wkb
  if (rkstage == 1)
    dxRay = 0.0
    dkRay = 0.0
    ddxRay = 0.0
  end


  cgx_max = 0.0
  cgy_max = 0.0
  cgz_max = 0.0
  cgz_max_tfc .= 0.0


  kz0 = ifelse(case_wkb == 3, 0, 1)

  for kz in kz0:nz, jy = 1:ny, ix = 1:nx
    ijk = CartesianIndex(ix, jy, kz)
    nskip = 0

    if (nray[ijk] < 1)
      continue
    end

    for iray in 1:nray(ix, jy, kz)
      rijk = CartesianIndex(iray, ijk)
      wnrk, wnrl, wnrm = wavenumbers(rijk, rays)
      # wnrk = rays.k[rijk]
      # wnrl = rays.l[rijk]
      # wnrm = rays.m[rijk]

      wnrh = sqrt(wnrk^2 + wnrl^2)

      xr, yr, zr = pos(rijk, rays)
      dxr, dry, dzr = extents(rijk, rays)

      # !skip ray volumes that have left the domain
      if (case_wkb != 3)
        if (zr1 < ztildetfc[ix, jy, -1])
          nskip = nskip + 1
          continue
        end
      end

      nnr = stratification(zr, 1 )
      nnr1 = stratification(zr1, 1)
      nnr2 = stratification(zr2, 1)

      omir1 =
        branchr * sqrt(nnr1 * wnrh^2 + f_cor_nd^2 * wnrm^2) /
        sqrt(wnrh^2 + wnrm^2)

      omir =
        branchr * sqrt(nnr * wnrh^2 + f_cor_nd^2 * wnrm^2) /
        sqrt(wnrh^2 + wnrm^2)

      if (nnr2 <= 0.0)
        # print *, 'NNr2 =', NNr2, '<= 0.0 at'
        # print *, 'zr2 =', zr2, 'from'
        # print *, 'ray(iRay,ix,jy,kz)%z =', rays[rijk]%z
        # print *, 'ray(iRay,ix,jy,kz)%dzray =', rays[rijk]%dzray
        # print *, 'iRay,ix,jy,kz =', rijk
        exit()
      end

      if (wnrh <= 0.0)
        # print *, 'wnrh =', wnrh, '<= 0.0 from'
        # print *, 'wnrk =', wnrk
        # print *, 'wnrl =', wnrl
        # print *, 'iRay,ix,jy,kz =', rijk
        exit()
      end

      omir2 =
        branchr * sqrt(nnr2 * wnrh^2 + f_cor_nd^2 * wnrm^2) /
        sqrt(wnrh^2 + wnrm^2)

      rays.omega[rijk] = omir

      # ! intrinsic group velocities at the respective edges of
      # ! the ray volumes

      if (sizex > 1)
        # ! intrinsic group velocity in x direction not depending
        # ! on x
        cgirx = wnrk * (nnr - omir^2) / (omir * (wnrh^2 + wnrm^2))
      end

      if (sizey > 1)
        # ! intrinsic group velocity in y direction not depending
        # ! on y
        cgiry = wnrl * (nnr - omir^2) / (omir * (wnrh^2 + wnrm^2))
      end

      # ! intrinsic vertical group velocity depending on z
      # ! (via the stratification)
      cgirz1 = -wnrm * (omir1^2 - f_cor_nd^2) / (omir1 * (wnrh^2 + wnrm^2))
      cgirz2 = -wnrm * (omir2^2 - f_cor_nd^2) / (omir2 * (wnrh^2 + wnrm^2))

      cgx_max, cgy_max, cgz_max = displacement(rijk, (cgx_max, cgy_max, cgz_max), cgz_max_tfc, rays, domain)

      # !-------------------------------
      # !    change of wavenumber
      # !-------------------------------

      # ! wave refraction only above lz(0) + zmin_wkb
      if (zr > lz[0] + zmin_wkb)
        #   # ! RK procedure

        # TODO: we updated xr etc. what does fortran do here? make a copy or reference?
        meanflow(xr, yr, zr, var, 11, dudxr)
        meanflow(xr, yr, zr, var, 12, dudyr)
        meanflow(xr, yr, zr, var, 13, dudzr)

        meanflow(xr, yr, zr, var, 21, dvdxr)
        meanflow(xr, yr, zr, var, 22, dvdyr)
        meanflow(xr, yr, zr, var, 23, dvdzr)

        if (zr < lz[0] - dz)
          # print *, 'ERROR IN setup_wkb: LOWER EDGE OF RAY  VOLUME', &
          #     &rijk, 'TOO LOW'
          exit()
        end

        if (zr < lz[0] - dz)
          # print *, 'ERROR IN transport_rayvol: RAY VOLUME', iRay, ix, &
          #     &jy, kz, 'TOO LOW'
          exit()
        end

        stratification(zr, 2, dnndzr)

        dkdt = -dudxr * wnrk - dvdxr * wnrl
        dldt = -dudyr * wnrk - dvdyr * wnrl
        dmdt =
          -dudzr * wnrk - dvdzr * wnrl -
          wnrh^2 * dnndzr / (2.0 * omir + (wnrh^2 + wnrm^2))

        dkRay[1, rijk] = dt * dkdt + alphaRK(rkStage) * dkRay[1, rijk]
        dkRay[2, rijk] = dt * dldt + alphaRK(rkStage) * dkRay[2, rijk]
        dkRay[3, rijk] = dt * dmdt + alphaRK(rkStage) * dkRay[3, rijk]

        rays.k[rijk] += betaRK(rkStage) * dkRay[1, rijk]
        rays.l[rijk] += betaRK(rkStage) * dkRay[2, rijk]
        rays.m[rijk] += betaRK[rkStage] * dkRay[3, rijk]

        # !----------------------------------------------
        # !    change of wave-number width of ray volumes
        # !----------------------------------------------

        # ! dk

        if (sizex > 1 && kz > 0 && !single_column)
          ddxdt = cgrx2 - cgrx1

          ddxRay[1, rijk] = dt * ddxdt + alphark[rkstage] * ddxray[1, rijk]

          rays.dxray[rijk] += betark[rkstage] * ddxray[1, rijk]

          if (rays.dxray[rijk] <= 0.0)
            # print *, 'dxray(', rijk, ') <= 0.0  ==> time &
            #     &step too large?'
            rays.dxray[rijk] *= -1
          end

          rays.dkray[rijk] = rays.area_xk[rijk] / rays.dxray[rijk]
        end

        # ! dl

        if (sizey > 1 && kz > 0 && !single_column)
          ddydt = cgry2 - cgry1

          ddxray[2, rijk] = dt * ddydt + alphark[rkstage] * ddxray[2, rijk]

          rays.dyray[rijk] += betark[rkstage] * ddxray[2, rijk]

          if (rays.dyray[rijk] <= 0.0)
            # print *, 'dyray(', rijk, ') <= 0.0  ==> time &
            #     &step too large?'
            rays.dyray[rijk] *= -1
          end

          rays.dlray[rijk] = rays.area_yl[rijk] / rays.dyray[rijk]
        end

        # !dm

        ddzdt = cgrz2 - cgrz1

        ddxRay[3, rijk] = dt * ddzdt + alphark[rkstage] * ddxray[3, rijk]

        rays.dzray[rijk] += betark[rkstage] * ddxray[3, rijk]

        if (rays.dzray[rijk] <= 0.0)
          # print *, 'dzray(', rijk, ') <= 0.0  ==> time step &
          #     &too large?'
          rays.dzray[rijk] *= -1
        end

        rays.dmray[rijk] = rays.area_zm[rijk] / rays.dzray[rijk]
      end

      # !-----------------------------------
      # ! update of the intrinsic frequency
      # !-----------------------------------

      wnrk = rays.k[rijk]
      wnrl = rays.l[rijk]
      wnrm = rays.m[rijk]

      wnrh = sqrt(wnrk^2 + wnrl^2)

      zr = rays.z[rijk]

      stratification(zr, 1, NNr)

      omir =
        branchr * sqrt(nnr * wnrh^2 + f_cor_nd^2 * wnrm^2) /
        sqrt(wnrh^2 + wnrm^2)

      rays.omega[rijk] = omir

    end # ray loop
    if (nskip > 0)
      # print *, nskip, 'r.v. skipped in transport_rayvol out of', &
      # &nRay(ix, jy, kz)
    end
  end # grid loop

  # ! Sponge layer
  if (spongelayer && unifiedsponge)
    for kz = 1:nz, jy = 1:ny, ix = 1:nx, iray = 1:nray(ix, jy, kz)
      rijk = CartesianIndex(iRay, ix, jy, kz)
      xr = ray.x[rijk]
      yr = ray.y[rijk]
      zr = ray.z[rijk]
      alphaSponge = 2.0 * interpolate_sponge(xr, yr, zr)
      betasponge = 1.0 / (1.0 + alphasponge * stepfrac(rkstage) * dt)
      rays.dens[rijk] *= betaSponge
    end
  end

  if (case_wkb == 3 && !steady_state)
    orographic_source(var, ray, time, stepfrac(rkstage) * dt)
  end
end


function transport_rayvol_steady()
  grid = state.grid
  for kz = 1:nz, jy = 1:ny, ix = 1:nx

    # ! Set ray-volume count.
    nray(ix, jy, kz) = nray(ix, jy, kz - 1)

    # ! Set up saturation computation.
    integral1 = 0.0
    integral2 = 0.0
    m2b2 = 0.0
    m2b2k2 = 0.0

    # ! Loop over ray volumes.
    for iRay = 1:nray(ix, jy, kz)

      # ! Prepare ray volume.
      ray[rijk] = ray[rijk-1]

      # ! Skip modes with zero wave-action density.
      if (ray(rijk - 1).dens == 0.0)
        continue
      end

      # ! Set vertical position (and extent).
      ray[rijk].z =
        zTildeTFC(ix, jy, kz - 1) +
        (ray(rijk - 1).z - zTildeTFC(ix, jy, kz - 2)) / jac(ix, jy, kz - 1) *
        jac(ix, jy, kz)
      ray[rijk].dzray =
        ray(rijk - 1).dzray * jac(ix, jy, kz) / jac(ix, jy, kz - 1)

      # ! Get horizontal wavenumbers.
      wnrk = rays[rijk].k
      wnrl = rays[rijk].l
      wnrh = sqrt(wnrk^2.0 + wnrl^2.0)

      # ! Set reference level.
      kz0 = max(1, kz - 1)

      # ! Compute vertical group velocity at the level below.
      stratification(ray(rijk0).z, 1, NN_nd)
      omir = ray(rijk0).omega
      if (branchr * omir > f_cor_nd && branchr * omir < sqrt(NN_nd))
        wnrm = ray(rijk0).m
        cgirz0 =
          wnrm * (f_cor_nd^2 - NN_nd) * wnrh^2 / omir / (wnrh^2 + wnrm^2)^2
      else
        ray(rijk0).dens = 0.0
        rays[rijk].dens = 0.0
        continue
      end

      # ! Compute local intrinsic frequency, vertical
      # ! wavenumber and vertical group velocity.
      stratification(rays[rijk].z, 1, NN_nd)
      omir =
        -0.5 * (var.u(ix, jy, kz) + var.u(ix - 1, jy, kz)) * wnrk -
        0.5 * (var.v(ix, jy, kz) + var.v(ix, jy - 1, kz)) * wnrl
      if (branchr * omir > f_cor_nd && branchr * omir < sqrt(NN_nd))
        wnrm =
          -branchr * sqrt(wnrh^2 * (NN_nd - omir^2) / (omir^2 - f_cor_nd^2))
        cgirz =
          wnrm * (f_cor_nd^2 - NN_nd) * wnrh^2 / omir / (wnrh^2 + wnrm^2)^2
      else
        ray(rijk0).dens = 0.0
        rays[rijk].dens = 0.0
        continue
      end

      # ! Set local intrinsic frequency and vertical wavenumber.
      rays[rijk].omega = omir
      rays[rijk].m = wnrm

      # ! Set local wave action density.
      if (spongeLayer && unifiedSponge)
        xr = rays[rijk].x
        yr = rays[rijk].y
        zr = rays[rijk].z
        alphaSponge = 2.0 * interpolate_sponge(xr, yr, zr)
        rays[rijk].dens =
          1.0 / (1.0 + alphaSponge / cgirz * (rays[rijk].z - ray(rijk0).z)) *
          cgirz0 *
          ray(rijk0).dens / cgirz
      else
        rays[rijk].dens = cgirz0 * ray(rijk0).dens / cgirz
      end

      # ! Cycle if saturation scheme is turned off.
      if (!lsaturation)
        continue
      end

      # ! Get ray volume extents.
      dxr = rays[rijk].dxray
      dyr = rays[rijk].dyray
      dzr = rays[rijk].dzray
      dwnrk = rays[rijk].dkray
      dwnrl = rays[rijk].dlray
      dwnrm = rays[rijk].dmray

      # ! Compute phase space factor.
      dzi = min(dzr, jac(ix, jy, kz) * dz)
      facpsp = dzi / jac(ix, jy, kz) / dz * dwnrm

      if (sizeX > 1)
        dxi = min(dxr, dx)
        facpsp = facpsp * dxi / dx * dwnrk
      end
      if (sizeY > 1)
        dyi = min(dyr, dy)
        facpsp = facpsp * dyi / dy * dwnrl
      end

      # ! Update saturation amplitude.
      integral1 = wnrh^2 * wnrm^2 / ((wnrh^2 + wnrm^2) * omir) * facpsp
      m2b2 =
        m2b2 +
        2.0 * NN_nd^2 / rhoStratTFC(ix, jy, kz) * integral1 * rays[rijk].dens

      integral2 = wnrh^2 * wnrm^2 / omir * facpsp
      m2b2k2 =
        m2b2k2 +
        2.0 * NN_nd^2 / rhoStratTFC(ix, jy, kz) *
        integral2 *
        rays[rijk].dens *
        jac(ix, jy, kz) *
        dz / cgirz

    end
    # ! Compute diffusion coefficient
    stratification(zTFC(ix, jy, kz), 1, NN_nd)
    if (m2b2k2 == 0.0 || m2b2 < alpha_sat^2 * NN_nd^2)
      diffusion = 0.0
    else
      diffusion = (m2b2 - alpha_sat^2 * NN_nd^2) / (2.0 * m2b2k2)
    end
    # ! Reduce wave action density.
    for iray = 1:nray(ix, jy, kz)
      if (!lsaturation)
        continue
      end
      if (rays[rijk].dens == 0.0)
        continue
      end
      wnrk = rays[rijk].k
      wnrl = rays[rijk].l
      wnrm = rays[rijk].m
      wnrh = sqrt(wnrk^2 + wnrl^2)
      stratification(rays[rijk].z, 1, NN_nd)
      omir = rays[rijk].omega
      if (branchr * omir > f_cor_nd && branchr * omir < sqrt(NN_nd))
        cgirz =
          wnrm * (f_cor_nd^2 - NN_nd) * wnrh^2 / omir / (wnrh^2 + wnrm^2)^2
      else
        ray(rijk0).dens = 0.0
        rays[rijk].dens = 0.0
        cycle
      end
      rays[rijk].dens =
        rays[rijk].dens * max(
          0.0,
          1.0 -
          jac(ix, jy, kz) * dz / cgirz * 2.0 * diffusion * (wnrh^2 + wnrm^2),
        )

    end
  end
  if (sizeX > 1)
    setboundary_rayvol_x(ray)
  end
  if (sizeY > 1)
    setboundary_rayvol_y(ray)
  end
  return
end


function transport_rayvol()
  # preamble
  #
  # ix0 = is + nbx - 1
  # jy0 = js + nby - 1

  # f_cor_nd = f_Coriolis_dim * tRef

  # if(case_wkb == 3 .and. steady_state) then
  #   call orographic_source(var, ray, time, stepFrac(RKStage) * dt)
  # end if


  # big if steady_state(), then return
  #

end
