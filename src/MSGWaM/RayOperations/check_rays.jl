"""
```julia
check_rays(state::State)
```

Check if all ray volumes are assigned to the correct grid cells.

# Arguments

  - `state`: Model state.
"""
function check_rays end

function check_rays(state::State)
    (; ndx, ndy) = state.namelists.domain
    (; io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, x, y, ztildetfc) = state.grid
    (; nray, rays) = state.wkb

    # Loop over ray volumes.
    @ivy for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

        for r in 1:nray[i, j, k]
            if rays.dens[r, i, j, k] == 0
                continue
            end

            # Check zonal position.
            if ndx > 1
                xr = rays.x[r, i, j, k]

                if xr < x[io + i] - dx / 2
                    println("Error in check_rays:")
                    println(
                        "xr = ",
                        xr,
                        " < x[io + i] - dx / 2 = ",
                        x[io + i] - dx / 2,
                    )
                    println("io = ", io)
                    println("(r, i, j, k) = ", (r, i, j, k))
                    exit()
                end

                if xr > x[io + i] + dx / 2
                    println("Error in check_rays:")
                    println(
                        "xr = ",
                        xr,
                        " > x[io + i] + dx / 2 = ",
                        x[io + i] + dx / 2,
                    )
                    println("io = ", io)
                    println("(r, i, j, k) = ", (r, i, j, k))
                    exit()
                end
            end

            # Check meridional position.
            if ndy > 1
                yr = rays.y[r, i, j, k]

                if yr < y[jo + j] - dy / 2
                    println("Error in check_rays:")
                    println(
                        "yr = ",
                        yr,
                        " < y[jo + j] - dy / 2 = ",
                        y[jo + j] - dy / 2,
                    )
                    println("jo = ", jo)
                    println("(r, i, j, k) = ", (r, i, j, k))
                    exit()
                end

                if yr > y[jo + j] + dy / 2
                    println("Error in check_rays:")
                    println(
                        "yr = ",
                        yr,
                        " > y[jo + j] + dy / 2 = ",
                        y[jo + j] + dy / 2,
                    )
                    println("jo = ", jo)
                    println("(r, i, j, k) = ", (r, i, j, k))
                    exit()
                end
            end

            # Check vertical position.
            zr = rays.z[r, i, j, k]

            if zr < ztildetfc[i, j, k - 1]
                println("Error in check_rays:")
                println(
                    "zr = ",
                    zr,
                    " < ztildetfc[i, j, k - 1] = ",
                    ztildetfc[i, j, k - 1],
                )
                println("(r, i, j, k) = ", (r, i, j, k))
                exit()
            end

            if zr > ztildetfc[i, j, k]
                println("Error in check_rays:")
                println(
                    "zr = ",
                    zr,
                    " > ztildetfc[i, j, k] = ",
                    ztildetfc[i, j, k],
                )
                println("(r, i, j, k) = ", (r, i, j, k))
                exit()
            end
        end
    end

    return
end
