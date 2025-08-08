"""
```julia
compute_horizontal_cell_indices(
    state::State,
    xr::AbstractFloat,
    yr::AbstractFloat,
    dxr::AbstractFloat,
    dyr::AbstractFloat,
)
```

From the given horizontal ray-volume position and extent, determine the indices of the grid cells that contain the ray-volume edges.

# Arguments

- `state`: Model state.
- `xr`: Ray-volume position in ``x``.
- `yr`: Ray-volume position in ``y``.
- `dxr`: Ray-volume extent in ``x``.
- `dyr`: Ray-volume extent in ``y``.

# Returns

- `::Integer`: Zonal index of the grid cell that contains the left edge.
- `::Integer`: Zonal index of the grid cell that contains the right edge.
- `::Integer`: Meridional index of the grid cell that contains the backward edge.
- `::Integer`: Meridional index of the grid cell that contains the forward edge.
"""
function compute_horizontal_cell_indices end

function compute_horizontal_cell_indices(
    state::State,
    xr::AbstractFloat,
    yr::AbstractFloat,
    dxr::AbstractFloat,
    dyr::AbstractFloat,
)
    (; sizex, sizey) = state.namelists.domain
    (; i0, i1, j0, j1, io, jo) = state.domain
    (; lx, ly, dx, dy) = state.grid

    if sizex > 1
        ixmin = floor(Int, (xr - lx[1] - dxr / 2) / dx) + i0 - io
        ixmax = floor(Int, (xr - lx[1] + dxr / 2) / dx) + i0 - io

        if ixmin > i1 + 1
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
        ixmin = i0
        ixmax = i0
    end

    if sizey > 1
        jymin = floor(Int, (yr - ly[1] - dyr / 2) / dy) + j0 - jo
        jymax = floor(Int, (yr - ly[1] + dyr / 2) / dy) + j0 - jo

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
        jymin = j0
        jymax = j0
    end

    return (ixmin, ixmax, jymin, jymax)
end
