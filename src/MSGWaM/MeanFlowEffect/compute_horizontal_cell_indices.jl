"""
    compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

Calculate grid cell indices that overlap with a ray volume in horizontal dimensions.

# Arguments

  - `state`: Simulation state containing domain and grid information
  - `xr`: Ray center x-coordinate in computational space
  - `yr`: Ray center y-coordinate in computational space
  - `dxr`: Ray volume width in x-direction
  - `dyr`: Ray volume width in y-direction

# Returns

  - `Tuple{Int,Int,Int,Int}`: (ixmin, ixmax, jymin, jymax) - minimum and maximum grid indices
    in x and y directions that overlap with the ray volume

# Implementation

For each direction, calculates cell indices by:

 1. Finding grid cells containing ray volume boundaries
 2. Ensuring indices remain within computational domain
 3. Handling 1D domains by setting min/max indices to the same value

# Error Checking

Throws error if ray volume falls completely outside computational domain
"""
function compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)
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
