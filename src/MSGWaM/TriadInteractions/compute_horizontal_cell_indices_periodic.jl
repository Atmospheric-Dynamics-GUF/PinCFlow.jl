function compute_horizontal_cell_indices_periodic end

function compute_horizontal_cell_indices_periodic(
    state::State,
    xr::AbstractFloat,
    yr::AbstractFloat,
    dxr::AbstractFloat,
    dyr::AbstractFloat,
)::NTuple{4, <:Integer}

    (; x_size, y_size) = state.namelists.domain
    (; i0, j0, io, jo, nxx, nyy) = state.domain
    (; lx, ly, dx, dy) = state.grid

    if x_size > 1
        imin = floor(Int, (xr + lx / 2 - dxr / 2) / dx) + i0 - io
        imax = floor(Int, (xr + lx / 2 + dxr / 2) / dx) + i0 - io

        if imin < 1 || imax > nxx
            error("Ray volume exceeds available x halo: imin = $imin, imax = $imax, valid range = 1:$nxx")
        end
    else
        imin = i0
        imax = i0
    end

    if y_size > 1
        jmin = floor(Int, (yr + ly / 2 - dyr / 2) / dy) + j0 - jo
        jmax = floor(Int, (yr + ly / 2 + dyr / 2) / dy) + j0 - jo

        if jmin < 1 || jmax > nyy
            error("Ray volume exceeds available y halo: jmin = $jmin, jmax = $jmax, valid range = 1:$nyy")
        end
    else
        jmin = j0
        jmax = j0
    end

    return (imin, imax, jmin, jmax)
end