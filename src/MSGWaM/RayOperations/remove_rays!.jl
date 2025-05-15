function remove_rays!(state::State)
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; nray, rays) = state.wkb

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for kz in kz0:kz1, jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        nrlc = 0
        for iray in 1:nray[ix, jy, kz]
            if rays.dens[iray, ix, jy, kz] == 0
                continue
            end
            nrlc += 1
            if nrlc != iray
                copy_rays!(rays, (iray, ix, jy, kz), (nrlc, ix, jy, kz))
                rays.dens[iray, ix, jy, kz] = 0.0
            end
        end
        nray[ix, jy, kz] = nrlc
    end

    return
end
