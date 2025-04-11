function calc_gwh!(state::State)
    (; gwh) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; ztfc) = state.grid
    (; sizex, sizey, lz) = state.namelists.domain
    (; integrals) = state.wkb
    (; f_cor_nd) = state.namelists.atmosphere

    gwh .= 0.0

    if f_cor_nd ≠ 0.0 && (sizex > 1 || sizey > 1)
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            if ztfc[ix, jy, kz] < (lz[1] + zmin_wkb)
                continue
            end

            gwh[ix, jy, kz] = integrals.dthetadt[ix, jy, kz]
        end
    end
    return
end
