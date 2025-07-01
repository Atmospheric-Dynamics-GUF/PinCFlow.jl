"""
```julia
check_rays(state::State)
```

Validate ray volume positions are within expected grid cell boundaries.

Performs consistency checks to ensure all ray volumes are positioned
within their assigned grid cells and haven't drifted outside acceptable
bounds during propagation or other operations.

# Arguments

  - `state::State`: Complete simulation state containing ray data

# Validation Checks

## Zonal Position (X-direction)

For multi-dimensional domains (`sizex > 1`):

  - Ray x-position must be within: `[x_cell - dx/2, x_cell + dx/2]`
  - Accounts for grid cell center and spacing

## Meridional Position (Y-direction)

For multi-dimensional domains (`sizey > 1`):

  - Ray y-position must be within: `[y_cell - dy/2, y_cell + dy/2]`
  - Validates against cell boundaries in y-direction

## Vertical Position (Z-direction)

For all ray volumes:

  - Must satisfy: `z_bottom ≤ ray_z ≤ z_top`
  - Uses terrain-following coordinate levels
  - Checks against `ztildetfc` level boundaries

# Grid Cell Range

Validates rays in extended domain including halo regions:

  - X: `(i0-1):(i1+1)`
  - Y: `(j0-1):(j1+1)`
  - Z: `(k0-1):(k1+1)`

# Error Handling

  - Prints detailed error information if violations found
  - Reports ray index, grid cell coordinates, and constraint values
  - Calls `exit()` to terminate simulation on validation failure

# Use Cases

  - **Debugging**: Identify ray propagation issues
  - **Code verification**: Ensure algorithms maintain ray locality
  - **Performance**: Detect excessive ray spreading    # Loop over ray volumes.

# Performance Impact

This is a validation function that can be expensive for large ray counts.
Typically used during development/debugging rather than production runs.
"""
function check_rays(state::State)
    (; sizex, sizey) = state.namelists.domain
    (; io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, x, y, ztildetfc) = state.grid
    (; nray, rays) = state.wkb

    # Loop over ray volumes.
    for kz in (k0 - 1):(k1 + 1),
        jy in (j0 - 1):(j1 + 1),
        ix in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ix, jy, kz]
            if rays.dens[iray, ix, jy, kz] == 0
                continue
            end

            # Check zonal position.
            if sizex > 1
                xr = rays.x[iray, ix, jy, kz]

                if xr < x[io + ix] - dx / 2
                    println("Error in check_rays:")
                    println(
                        "xr = ",
                        xr,
                        " < x[io + ix] - dx / 2 = ",
                        x[io + ix] - dx / 2,
                    )
                    println("io = ", io)
                    println("(iray, ix, jy, kz) = ", (iray, ix, jy, kz))
                    exit()
                end

                if xr > x[io + ix] + dx / 2
                    println("Error in check_rays:")
                    println(
                        "xr = ",
                        xr,
                        " > x[io + ix] + dx / 2 = ",
                        x[io + ix] + dx / 2,
                    )
                    println("io = ", io)
                    println("(iray, ix, jy, kz) = ", (iray, ix, jy, kz))
                    exit()
                end
            end

            # Check meridional position.
            if sizey > 1
                yr = rays.y[iray, ix, jy, kz]

                if yr < y[jo + jy] - dy / 2
                    println("Error in check_rays:")
                    println(
                        "yr = ",
                        yr,
                        " < y[jo + jy] - dy / 2 = ",
                        y[jo + jy] - dy / 2,
                    )
                    println("jo = ", jo)
                    println("(iray, ix, jy, kz) = ", (iray, ix, jy, kz))
                    exit()
                end

                if yr > y[jo + jy] + dy / 2
                    println("Error in check_rays:")
                    println(
                        "yr = ",
                        yr,
                        " > y[jo + jy] + dy / 2 = ",
                        y[jo + jy] + dy / 2,
                    )
                    println("jo = ", jo)
                    println("(iray, ix, jy, kz) = ", (iray, ix, jy, kz))
                    exit()
                end
            end

            # Check vertical position.
            zr = rays.z[iray, ix, jy, kz]

            if zr < ztildetfc[ix, jy, kz - 1]
                println("Error in check_rays:")
                println(
                    "zr = ",
                    zr,
                    " < ztildetfc[ix, jy, kz - 1] = ",
                    ztildetfc[ix, jy, kz - 1],
                )
                println("(iray, ix, jy, kz) = ", (iray, ix, jy, kz))
                exit()
            end

            if zr > ztildetfc[ix, jy, kz]
                println("Error in check_rays:")
                println(
                    "zr = ",
                    zr,
                    " > ztildetfc[ix, jy, kz] = ",
                    ztildetfc[ix, jy, kz],
                )
                println("(iray, ix, jy, kz) = ", (iray, ix, jy, kz))
                exit()
            end
        end
    end

    return
end
