"""
```julia
compute_horizontal_cell_indices(
    state::State,
    xr::AbstractFloat,
    yr::AbstractFloat,
    dxr::AbstractFloat,
    dyr::AbstractFloat,
)::NTuple{4, <:Integer}
```

From the given horizontal ray-volume position and extent, determine the indices of the grid cells that contain the ray-volume edges and return them (in the order left, right, backward and forward).

# Arguments

  - `state`: Model state.

  - `xr`: Ray-volume position in ``x``.

  - `yr`: Ray-volume position in ``y``.

  - `dxr`: Ray-volume extent in ``x``.

  - `dyr`: Ray-volume extent in ``y``.
"""
function compute_horizontal_cell_indices end

function compute_horizontal_cell_indices(
    state::State,
    xr::AbstractFloat,
    yr::AbstractFloat,
    dxr::AbstractFloat,
    dyr::AbstractFloat,
)::NTuple{4, <:Integer}
    (; x_size, y_size) = state.namelists.domain
    (; i0, i1, j0, j1, io, jo) = state.domain
    (; lx, ly, dx, dy) = state.grid

    if x_size > 1
        imin = floor(Int, (xr + lx / 2 - dxr / 2) / dx) + i0 - io
        imax = floor(Int, (xr + lx / 2 + dxr / 2) / dx) + i0 - io

        if imin > i1 + 1
            error(
                "Error in meanflow_effect: imin = ",
                imin,
                " > i1 + 1 = ",
                i1 + 1,
            )
        else
            imin = max(imin, i0)
        end

        if imax < i0 - 1
            error(
                "Error in meanflow_effec: imax =",
                imax,
                "< i0 - 1 = ",
                i0 - 1,
            )
        else
            imax = min(imax, i1)
        end

    else
        imin = i0
        imax = i0
    end

    if y_size > 1
        jmin = floor(Int, (yr + ly / 2 - dyr / 2) / dy) + j0 - jo
        jmax = floor(Int, (yr + ly / 2 + dyr / 2) / dy) + j0 - jo

        if jmin > j1 + 1
            error(
                "Error in meanflow_effect: jmin = ",
                jmin,
                " > j1 + 1 = ",
                j1 + 1,
            )
        else
            jmin = max(jmin, j0)
        end

        if jmax < j0 - 1
            error(
                "Error in meanflow_effec: jmax =",
                jmax,
                "< j0 - 1 = ",
                j0 - 1,
            )
        else
            jmax = min(jmax, j1)
        end
    else
        jmin = j0
        jmax = j0
    end

    return (imin, imax, jmin, jmax)
end
