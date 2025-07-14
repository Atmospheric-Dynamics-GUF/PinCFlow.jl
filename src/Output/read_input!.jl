"""
```julia
read_input!(state::State)
```

Read initial values for all prognostic variables from an HDF5 input file.

# Arguments

  - `state`: Model state.
"""
function read_input!(state::State)

    # Get all necessary fields.
    (; sizex, sizey) = state.namelists.domain
    (; iin, input_file) = state.namelists.output
    (; model, testcase) = state.namelists.setting
    (; comm, sizezz, nx, ny, nz, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) =
        state.domain
    (; lref, tref, rhoref, uref) = state.constants
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

    # Open the file. Note: Fused in-place assignments cannot be used here!
    time = h5open(input_file, "r", comm) do file

        # Read the time.
        time = file["t"][iin] / tref

        # Read the density fluctuations.
        @views rhop[i0:i1, j0:j1, k0:k1] =
            file["rhop"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iin,
            ] ./ rhoref
        if model != Boussinesq()
            @views rho[i0:i1, j0:j1, k0:k1] .= rhop[i0:i1, j0:j1, k0:k1]
        end

        # Read the staggered zonal wind.
        @views u[i0:i1, j0:j1, k0:k1] =
            file["us"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iin,
            ] ./ uref

        # Read the staggered meridional wind.
        @views v[i0:i1, j0:j1, k0:k1] =
            file["vs"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iin,
            ] ./ uref

        # Read the staggered transformed vertical wind.
        @views w[i0:i1, j0:j1, k0:k1] =
            file["wstfc"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iin,
            ] ./ uref

        # Read the Exner-pressure fluctuations.
        @views pip[i0:i1, j0:j1, k0:k1] = file["pip"][
            (io + 1):(io + nx),
            (jo + 1):(jo + ny),
            (ko + 1):(ko + nz),
            iin,
        ]

        # Read the mass-weighted potential temperature.
        if model == Compressible()
            @views p[i0:i1, j0:j1, k0:k1] = file["p"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iin,
            ]
        end

        # Read ray-volume properties.
        if typeof(testcase) <: AbstractWKBTestCase
            dk0 = ko == 0 ? 1 : 0

            for (output_name, field_name) in zip(
                ("xr", "yr", "zr", "dxr", "dyr", "dzr"),
                (:x, :y, :z, :dxray, :dyray, :dzray),
            )
                @views getfield(rays, field_name)[
                    1:nray_max,
                    i0:i1,
                    j0:j1,
                    (k0 - dk0):k1,
                ] =
                    file[output_name][
                        1:nray_max,
                        (io + 1):(io + nx),
                        (jo + 1):(jo + ny),
                        (ko + 2 - dk0):(ko + nz + 1),
                        iin,
                    ] ./ lref
            end

            for (output_name, field_name) in zip(
                ("kr", "lr", "mr", "dkr", "dlr", "dmr"),
                (:k, :l, :m, :dkray, :dlray, :dmray),
            )
                @views getfield(rays, field_name)[
                    1:nray_max,
                    i0:i1,
                    j0:j1,
                    (k0 - dk0):k1,
                ] =
                    file[output_name][
                        1:nray_max,
                        (io + 1):(io + nx),
                        (jo + 1):(jo + ny),
                        (ko + 2 - dk0):(ko + nz + 1),
                        iin,
                    ] .* lref
            end

            @views rays.dens[1:nray_max, i0:i1, j0:j1, (k0 - dk0):k1] =
                file["nr"][
                    1:nray_max,
                    (io + 1):(io + nx),
                    (jo + 1):(jo + ny),
                    (ko + 2 - dk0):(ko + nz + 1),
                    iin,
                ] ./ rhoref ./ uref .^ 2 ./ tref ./ lref .^ dim

            # Determine nray.
            for kz in (k0 - dk0):k1, jy in j0:j1, ix in i0:i1
                nrlc = 0
                for iray in 1:nray_max
                    if rays.dens[iray, ix, jy, kz] == 0
                        continue
                    end
                    nrlc += 1
                end
                nray[ix, jy, kz] = nrlc
            end
        end

        # Return.
        return time
    end

    # Return.
    return time
end
