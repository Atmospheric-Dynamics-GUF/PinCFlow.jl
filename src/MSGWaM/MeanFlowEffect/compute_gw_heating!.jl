function compute_gw_heating!(state::State)
    (; sizex, sizey) = state.namelists.domain
    (; f_coriolis_dim) = state.namelists.atmosphere
    (; zmin_wkb_dim) = state.namelists.wkb
    (; lref, tref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; lz, ztfc) = state.grid
    (; gwh) = state.variables.predictands
    (; integrals) = state.wkb

    f_cor_nd = f_coriolis_dim * tref

    gwh .= 0.0

    if f_cor_nd != 0.0 && (sizex > 1 || sizey > 1)
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            if ztfc[ix, jy, kz] < (lz[1] + zmin_wkb_dim / lref)
                continue
            end

            gwh[ix, jy, kz] = integrals.dthetadt[ix, jy, kz]
        end
    end
    return
end
