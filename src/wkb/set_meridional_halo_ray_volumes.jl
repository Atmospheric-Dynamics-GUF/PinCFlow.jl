function set_meridional_halo_ray_volumes!(state::State)
  (; i0, i1, j0, j1, k0, k1, back, forw) = state.domain
  (; fields, send_rays_back, send_rays_forw, recv_rays_back, recv_rays_forw) =
    state.wkb.ray_communication
  (; nray, rays) = state.wkb

  # Loop over grid cells in z and x.
  for kz in (k0 - 1):(k1 + 1), ix in (i0 - 1):(i1 + 1)

    # Set backward boundary ray volumes.
    if nray[ix, j0 - 1, kz] > 0
      for iray in 1:nray[ix, j0 - 1, kz]
        for (index, field) in enumerate(fields)
          @views send_rays_forw[index] = getfield(rays, field)[iray, ix, j1, kz]
        end

        MPI.Sendrecv!(
          send_rays_forw,
          recv_rays_back,
          comm;
          dest = forw,
          source = back,
        )

        for (index, field) in enumerate(fields)
          @views getfield(rays, field)[iray, ix, j0 - 1, kz] =
            recv_rays_back[index]
        end
      end
    end

    # Set surfaces of backward boundary ray volumes.
    if nray[ix, j0 - 1, kz] > 0
      for iray in 1:nray[ix, j0 - 1, kz]
        rays.area_xk[iray, ix, j0 - 1, kz] =
          rays.dxray[iray, ix, j0 - 1, kz] * rays.dkray[iray, ix, j0 - 1, kz]
        rays.area_yl[iray, ix, j0 - 1, kz] =
          rays.dyray[iray, ix, j0 - 1, kz] * rays.dlray[iray, ix, j0 - 1, kz]
        rays.area_zm[iray, ix, j0 - 1, kz] =
          rays.dzray[iray, ix, j0 - 1, kz] * rays.dmray[iray, ix, j0 - 1, kz]
      end
    end

    # Set forward boundary ray volumes.
    if nray[ix, j1 + 1, kz] > 0
      for iray in 1:nray[ix, j1 + 1, kz]
        for (index, field) in enumerate(fields)
          @views send_rays_back[index] = getfield(rays, field)[iray, ix, j0, kz]
        end

        MPI.Sendrecv!(
          send_rays_back,
          recv_rays_forw,
          comm;
          dest = back,
          source = forw,
        )

        for (index, field) in enumerate(fields)
          @views getfield(rays, field)[iray, ix, j1 + 1, kz] =
            recv_rays_forw[index]
        end
      end
    end

    # Set surfaces of forward boundary ray volumes.
    if nray[ix, j1 + 1, kz] > 0
      for iray in 1:nray[ix, j1 + 1, kz]
        rays.area_xk[iray, ix, j1 + 1, kz] =
          rays.dxray[iray, ix, j1 + 1, kz] * rays.dkray[iray, ix, j1 + 1, kz]
        rays.area_yl[iray, ix, j1 + 1, kz] =
          rays.dyray[iray, ix, j1 + 1, kz] * rays.dlray[iray, ix, j1 + 1, kz]
        rays.area_zm[iray, ix, j1 + 1, kz] =
          rays.dzray[iray, ix, j1 + 1, kz] * rays.dmray[iray, ix, j1 + 1, kz]
      end
    end
  end

  return
end
