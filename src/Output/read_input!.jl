"""
```julia
read_input!(state::State)
```

Read initial values for all prognostic variables from an HDF5 input file.

# Arguments

  - `state`: Model state.
"""
function read_input! end

function read_input!(state::State)

    # Get all necessary fields.
    (; sizex, sizey) = state.namelists.domain
    (; iin, input_file) = state.namelists.output
    (; model, testcase) = state.namelists.setting
    (; comm, nx, ny, nz, io, jo, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; lref, tref, rhoref, uref, thetaref) = state.constants
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands
    (; nray_max, nray, rays) = state.wkb

    # Determine dimensionality.
    dim = 1
    if sizex > 1
        dim += 1
    end
    if sizey > 1
        dim += 1
    end

    # Define slices.
    dk0 = ko == 0 ? 1 : 0
    (r, i, j, k, kr) = (1:nray_max, i0:i1, j0:j1, k0:k1, (k0 - dk0):k1)
    (id, jd, kd, krd) = (
        (io + 1):(io + nx),
        (jo + 1):(jo + ny),
        (ko + 1):(ko + nz),
        (ko + 2 - dk0):(ko + nz + 1),
    )

    # Open the file. Note: Fused in-place assignments cannot be used here!
    time = h5open(input_file, "r", comm) do file

        # Read the time.
        time = file["t"][iin] / tref

        # Read the density fluctuations.
        rhop[i, j, k] = file["rhop"][id, jd, kd, iin] ./ rhoref
        if model != Boussinesq()
            for kk in k, jj in j, ii in i
                rho[ii, jj, kk] = rhop[ii, jj, kk]
            end
        end

        # Read the staggered zonal wind.
        u[i, j, k] = file["us"][id, jd, kd, iin] ./ uref

        # Read the staggered meridional wind.
        v[i, j, k] = file["vs"][id, jd, kd, iin] ./ uref

        # Read the staggered transformed vertical wind.
        w[i, j, k] = file["wstfc"][id, jd, kd, iin] ./ uref

        # Read the Exner-pressure fluctuations.
        pip[i, j, k] = file["pip"][id, jd, kd, iin]

        # Read the mass-weighted potential temperature.
        if model == Compressible()
            p[i, j, k] = file["p"][id, jd, kd, iin] ./ rhoref ./ thetaref
        end

        # Read ray-volume properties.
        if typeof(testcase) <: AbstractWKBTestCase
            for (output_name, field_name) in zip(
                ("xr", "yr", "zr", "dxr", "dyr", "dzr"),
                (:x, :y, :z, :dxray, :dyray, :dzray),
            )
                getfield(rays, field_name)[r, i, j, kr] =
                    file[output_name][r, id, jd, krd, iin] ./ lref
            end

            for (output_name, field_name) in zip(
                ("kr", "lr", "mr", "dkr", "dlr", "dmr"),
                (:k, :l, :m, :dkray, :dlray, :dmray),
            )
                getfield(rays, field_name)[r, i, j, kr] =
                    file[output_name][r, id, jd, krd, iin] .* lref
            end

            rays.dens[r, i, j, kr] =
                file["nr"][r, id, jd, krd, iin] ./ rhoref ./ uref .^ 2 ./ tref ./ lref .^ dim

            # Determine nray.
            for kk in kr, jj in j, ii in i
                nrlc = 0
                for rr in r
                    if rays.dens[rr, ii, jj, kk] == 0
                        continue
                    end
                    nrlc += 1
                end
                nray[ii, jj, kk] = nrlc
            end
        end

        # Return.
        return time
    end

    # Return.
    return time
end
