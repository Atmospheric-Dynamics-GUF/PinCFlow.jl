function read_input!(state::State)

    # Get all necessary fields.
    (; sizex, sizey) = state.namelists.domain
    (; iin, folder) = state.namelists.output
    (; testcase) = state.namelists.setting
    (; comm, nx, ny, nz, io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; lref, tref, rhoref, uref) = state.constants
    (; rho, rhop, u, v, w, pip) = state.variables.predictands
    (; rays, nray_max) = state.wkb

    # Determine dimensionality.
    dim = 1
    if sizex > 1
        dim += 1
    end
    if sizey > 1
        dim += 1
    end

    # Open the file. Note: Fused in-place assignments cannot be used here!
    time = h5open(folder * "/pincflow_input.h5", "r", comm) do file

        # Read the time.
        time = file["t"][iin] / tref

        # Read the density fluctuations.
        @views rhop[i0:i1, j0:j1, k0:k1] =
            file["rhop"][(io + 1):(io + nx), (jo + 1):(jo + ny), 1:nz, iin] ./ rhoref
        @views rho[i0:i1, j0:j1, k0:k1] .= rhop[i0:i1, j0:j1, k0:k1]

        # Read the staggered zonal wind.
        @views u[i0:i1, j0:j1, k0:k1] =
            file["us"][(io + 1):(io + nx), (jo + 1):(jo + ny), 1:nz, iin] ./
            uref

        # Read the staggered meridional wind.
        @views v[i0:i1, j0:j1, k0:k1] =
            file["vs"][(io + 1):(io + nx), (jo + 1):(jo + ny), 1:nz, iin] ./
            uref

        # Read the staggered transformed vertical wind.
        @views w[i0:i1, j0:j1, k0:k1] =
            file["wstfc"][(io + 1):(io + nx), (jo + 1):(jo + ny), 1:nz, iin] ./ uref

        # Read the Exner-pressure fluctuations.
        @views pip[i0:i1, j0:j1, k0:k1] =
            file["pip"][(io + 1):(io + nx), (jo + 1):(jo + ny), 1:nz, iin]

        # Read ray-volume properties.
        if typeof(testcase) <: AbstractWKBTestCase
            for (output_name, field_name) in zip(
                ("xr", "yr", "zr", "dxr", "dyr", "dzr"),
                (:x, :y, :z, :dxray, :dyray, :dzray),
            )
                @views getfield(rays, field_name)[1:nray_max, i0:i1, j0:j1, k0:k1] =
                    file[output_name][
                        1:nray_max,
                        (io + 1):(io + nx),
                        (jo + 1):(jo + ny),
                        1:nz,
                        iin,
                    ] ./ lref
            end

            for (output_name, field_name) in zip(
                ("kr", "lr", "mr", "dkr", "dlr", "dmr"),
                (:k, :l, :m, :dkray, :dlray, :dmray),
            )
                @views getfield(rays, field_name)[1:nray_max, i0:i1, j0:j1, k0:k1] =
                    file[output_name][
                        1:nray_max,
                        (io + 1):(io + nx),
                        (jo + 1):(jo + ny),
                        1:nz,
                        iin,
                    ] .* lref
            end

            @views rays.dens[1:nray_max, i0:i1, j0:j1, k0:k1] =
                file["nr"][
                    1:nray_max,
                    (io + 1):(io + nx),
                    (jo + 1):(jo + ny),
                    1:nz,
                    iin,
                ] ./ rhoref ./ uref .^ 2 ./ tref ./ lref .^ dim
        end

        # Return.
        return time
    end

    # Return.
    return time
end
