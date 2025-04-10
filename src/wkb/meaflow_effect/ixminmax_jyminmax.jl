
function ixminmax_jyminmax(state, xr, yr, dxr, dyr)
    (; i0, i1, j0, j1, io, jo) = state.domain
    (; lx, ly, dx, dy) = state.grid
    (; sizex, sizey) = state.namelist.domain

    # indices of range of cells touched by a ray volume 
    if sizex > 1
        ixmin = floor((xr - dxr * 0.5 - lx[1]) / dx) + i0 - io
        ixmax = floor((xr + dxr * 0.5 - lx[1]) / dx) + i0 - io

        if ixmin > (i1 + 1)
            error(
                "Error in meanflow_effect: ixmin = ",
                ixmin,
                " > i1 + 1 = ",
                i1 + 1,
            )
        else
            ixmin = max(ixmin, i0)
        end

        if ixmax < i0 - 1
            error(
                "Error in meanflow_effec: ixmax =",
                ixmax,
                "< i0 - 1 = ",
                i0 - 1,
            )
        else
            ixmax = min(ixmax, i1)
        end

    else
        ixmin = 1
        ixmax = 1
    end

    if sizey > 1
        jymin = floor((yr - dyr * 0.5 - ly[1]) / dy) + j0 - jo
        jymax = floor((yr + dyr * 0.5 - ly[1]) / dy) + j0 - jo

        if jymin > j1 + 1
            error(
                "Error in meanflow_effect: jymin = ",
                jymin,
                " > j1 + 1 = ",
                j1 + 1,
            )
        else
            jymin = max(jymin, j0)
        end

        if jymax < j0 - 1
            error(
                "Error in meanflow_effec: jymax =",
                jymax,
                "< j0 - 1 = ",
                j0 - 1,
            )
        else
            jymax = min(jymax, j1)
        end
    else
        jymin = 1
        jymax = 1
    end

    return ixmin, ixmax, jymin, jaymax
end
