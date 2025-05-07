function set_vertical_boundary_rays!(
    state::State,
    zboundaries::SolidWallBoundaries,
)
    (; namelists, domain) = state
    (; npz) = namelists.domain
    (; sizezz, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; lx, ly, lz, dx, dy, topography_surface) = state.grid
    (; nray, rays) = state.wkb

    # Set ray-volume count and ray-volumes properties.
    if npz > 1
        set_vertical_halos_of_field!(
            nray,
            namelists,
            domain,
            zboundaries;
            layers = (1, 1, 1),
        )
        set_vertical_halo_rays!(state)
    end

    # Reflect ray volumes at the lower boundary.
    if ko == 0
        for kz in k0:(k0 + 1), jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
            for iray in 1:nray[ix, jy, kz]
                xr = rays.x[iray, ix, jy, kz]
                yr = rays.y[iray, ix, jy, kz]
                zr = rays.z[iray, ix, jy, kz]
                dzr = rays.dzray[iray, ix, jy, kz]
                wnrm = rays.m[iray, ix, jy, kz]

                ixrv = floor(Int, (xr - lx[1]) / dx) + i0 - io
                jyrv = floor(Int, (yr - ly[1]) / dy) + j0 - jo
                if topography_surface[ixrv, jyrv] - zr + 0.5 * dzr > eps()
                    rays.z[iray, ix, jy, kz] =
                        2.0 * topography_surface[ixrv, jyrv] - zr + dzr
                    rays.m[iray, ix, jy, kz] = -wnrm
                end
            end
        end
    end

    # Cut ray volumes at the upper boundary.
    if ko + nzz == sizezz
        for kz in (k1 - 1):k1, jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
            nrlc = 0
            for iray in 1:nray[ix, jy, kz]
                zr = rays.z[iray, ix, jy, kz]
                dzr = rays.dzray[iray, ix, jy, kz]

                if zr - 0.5 * dzr > lz[2]
                    continue
                end
                if zr + 0.5 * dzr > lz[2]
                    rays.dzray[iray, ix, jy, kz] = lz[2] - zr + 0.5 * dzr
                    rays.z[iray, ix, jy, kz] =
                        lz[2] - 0.5 * rays.dzray[iray, ix, jy, kz]
                end

                nrlc += 1
                if nrlc != iray
                    copy_rays!(rays, (iray, ix, jy, kz), (nrlc, ix, jy, kz))
                end
            end
            nray[ix, jy, kz] = nrlc
        end
    end

    return
end
