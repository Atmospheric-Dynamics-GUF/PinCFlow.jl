function set_zonal_halo_ray_volumes!(state::State)
  (; i0, i1, j0, j1, k0, k1, left, right) = state.domain
  (; fields, send_rays_left, send_rays_right, recv_rays_left, recv_rays_right) =
    state.wkb.ray_communication
  (; nray, rays) = state.wkb

  # Loop over grid cells in z and y.
  for kz in (k0 - 1):(k1 + 1), jy in (j0 - 1):(j1 + 1)

    # Set left boundary ray volumes.
    if nray[i0 - 1, jy, kz] > 0
      for iray in 1:nray[i0 - 1, jy, kz]
        for (index, field) in enumerate(fields)
          @views send_rays_right[index] =
            getfield(rays, field)[iray, i1, jy, kz]
        end

        MPI.Sendrecv!(
          send_rays_right,
          recv_rays_left,
          comm;
          dest = right,
          source = left,
        )

        for (index, field) in enumerate(fields)
          @views getfield(rays, field)[iray, i0 - 1, jy, kz] =
            recv_rays_left[index]
        end
      end
    end

    # Set surfaces of left boundary ray volumes.
    if nray[i0 - 1, jy, kz] > 0
      for iray in 1:nray[i0 - 1, jy, kz]
        rays.area_xk[iray, i0 - 1, jy, kz] =
          rays.dxray[iray, i0 - 1, jy, kz] * rays.dkray[iray, i0 - 1, jy, kz]
        rays.area_yl[iray, i0 - 1, jy, kz] =
          rays.dyray[iray, i0 - 1, jy, kz] * rays.dlray[iray, i0 - 1, jy, kz]
        rays.area_zm[iray, i0 - 1, jy, kz] =
          rays.dzray[iray, i0 - 1, jy, kz] * rays.dmray[iray, i0 - 1, jy, kz]
      end
    end

    # Set right boundary ray volumes.
    if nray[i1 + 1, jy, kz] > 0
      for iray in 1:nray[i1 + 1, jy, kz]
        for (index, field) in enumerate(fields)
          @views send_rays_left[index] = getfield(rays, field)[iray, i0, jy, kz]
        end

        MPI.Sendrecv!(
          send_rays_left,
          recv_rays_right,
          comm;
          dest = left,
          source = right,
        )

        for (index, field) in enumerate(fields)
          @views getfield(rays, field)[iray, i1 + 1, jy, kz] =
            recv_rays_right[index]
        end
      end
    end

    # Set surfaces of right boundary ray volumes.
    if nray[i1 + 1, jy, kz] > 0
      for iray in 1:nray[i1 + 1, jy, kz]
        rays.area_xk[iray, i1 + 1, jy, kz] =
          rays.dxray[iray, i1 + 1, jy, kz] * rays.dkray[iray, i1 + 1, jy, kz]
        rays.area_yl[iray, i1 + 1, jy, kz] =
          rays.dyray[iray, i1 + 1, jy, kz] * rays.dlray[iray, i1 + 1, jy, kz]
        rays.area_zm[iray, i1 + 1, jy, kz] =
          rays.dzray[iray, i1 + 1, jy, kz] * rays.dmray[iray, i1 + 1, jy, kz]
      end
    end
  end

  return
end
