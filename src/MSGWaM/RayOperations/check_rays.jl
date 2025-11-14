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
    (; x_size, y_size) = state.namelists.domain
    (; io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, x, y, zctilde) = state.grid
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
            if x_size > 1
                xr = rays.x[r, i, j, k]

                if xr < x[i] - dx / 2
                    error(
                        "Error in check_rays:\nxr = ",
                        xr,
                        " < x[i] - dx / 2 = ",
                        x[i] - dx / 2,
                        "\n(io, jo, ko) = ",
                        (io, jo, ko),
                        "\n(r, i, j, k) = ",
                        (r, i, j, k),
                    )
                end

                if xr > x[i] + dx / 2
                    error(
                        "Error in check_rays:\nxr = ",
                        xr,
                        " > x[i] + dx / 2 = ",
                        x[i] + dx / 2,
                        "\n(io, jo, ko) = ",
                        (io, jo, ko),
                        "\n(r, i, j, k) = ",
                        (r, i, j, k),
                    )
                end
            end

            # Check meridional position.
            if y_size > 1
                yr = rays.y[r, i, j, k]

                if yr < y[j] - dy / 2
                    error(
                        "Error in check_rays:\nyr = ",
                        yr,
                        " < y[j] - dy / 2 = ",
                        y[j] - dy / 2,
                        "\n(io, jo, ko) = ",
                        (io, jo, ko),
                        "\n(r, i, j, k) = ",
                        (r, i, j, k),
                    )
                end

                if yr > y[j] + dy / 2
                    error(
                        "Error in check_rays:\nyr = ",
                        yr,
                        " > y[j] + dy / 2 = ",
                        y[j] + dy / 2,
                        "\n(io, jo, ko) = ",
                        (io, jo, ko),
                        "\n(r, i, j, k) = ",
                        (r, i, j, k),
                    )
                end
            end

            # Check vertical position.
            zr = rays.z[r, i, j, k]

            if zr < zctilde[i, j, k - 1]
                error(
                    "Error in check_rays:\nzr = ",
                    zr,
                    " < zctilde[i, j, k - 1] = ",
                    zctilde[i, j, k - 1],
                    "\n(io, jo, ko) = ",
                    (io, jo, ko),
                    "\n(r, i, j, k) = ",
                    (r, i, j, k),
                )
            end

            if zr > zctilde[i, j, k]
                error(
                    "Error in check_rays:\nzr = ",
                    zr,
                    " > zctilde[i, j, k] = ",
                    zctilde[i, j, k],
                    "\n(io, jo, ko) = ",
                    (io, jo, ko),
                    "\n(r, i, j, k) = ",
                    (r, i, j, k),
                )
            end
        end
    end

    return
end
