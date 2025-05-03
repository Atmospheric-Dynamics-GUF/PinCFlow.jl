function write_output(
    state::State,
    time::AbstractFloat,
    iout::Integer,
    machine_start_time::DateTime,
)

    # Get all necessary fields.
    (; domain, grid) = state
    (; sizex, sizey, sizez) = state.namelists.domain
    (; prepare_restart, save_ray_volumes, output_variables, output_file) =
        state.namelists.output
    (; model, testcase) = state.namelists.setting
    (;
        comm,
        master,
        sizezz,
        nx,
        ny,
        nz,
        nzz,
        io,
        jo,
        ko,
        i0,
        i1,
        j0,
        j1,
        k0,
        k1,
    ) = domain
    (; tref, lref, rhoref, thetaref, uref) = state.constants
    (; x, y, ztfc) = grid
    (; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = state.atmosphere
    (; predictands) = state.variables
    (; rho, rhop, u, v, w, pip, p) = predictands
    (; nray_max, rays, tendencies) = state.wkb

    # Print information.
    if master
        println(repeat("-", 80))
        println("Output into file pincflow_output.h5")
        println("Physical time: ", time * tref, " s")
        println("Machine time: ", canonicalize(now() - machine_start_time))
        println(repeat("-", 80))
        println("")
    end

    # Advance output counter.
    iout += 1

    # Determine dimensionality.
    dim = 1
    if sizex > 1
        dim += 1
    end
    if sizey > 1
        dim += 1
    end

    # Open the file. Note: Fused in-place assignments cannot be used here!
    h5open(output_file, "r+", comm) do file

        # Write the time.
        HDF5.set_extent_dims(file["t"], (iout,))
        file["t"][iout] = time * tref

        # Write the horizontal grid.
        if iout == 1
            @views file["x"][:] = x[i0:(i0 + sizex - 1)] .* lref
            @views file["y"][:] = y[j0:(j0 + sizey - 1)] .* lref
        end

        # Write the vertical grid.
        if iout == 1
            @views file["z"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
            ] = ztfc[i0:i1, j0:j1, k0:k1] .* lref
        end

        # Write the background density.
        if model != Boussinesq() && iout == 1
            @views file["rhobar"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
            ] = rhostrattfc[i0:i1, j0:j1, k0:k1] .* rhoref
        end

        # Write the background potential temperature.
        if model != Boussinesq() && iout == 1
            @views file["thetabar"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
            ] = thetastrattfc[i0:i1, j0:j1, k0:k1] .* thetaref
        end

        # Write the buoyancy frequency.
        if model == Compressible()
            HDF5.set_extent_dims(file["n2"], (sizex, sizey, sizez, iout))
            @views file["n2"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] = bvsstrattfc[i0:i1, j0:j1, k0:k1] ./ tref .^ 2
        elseif model != Boussinesq() && iout == 1
            @views file["n2"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
            ] = bvsstrattfc[i0:i1, j0:j1, k0:k1] ./ tref .^ 2
        end

        # Write the mass-weighted potential temperature.
        if model == Compressible()
            HDF5.set_extent_dims(file["p"], (sizex, sizey, sizez, iout))
            @views file["p"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] = p[i0:i1, j0:j1, k0:k1] .* rhoref .* thetaref
        elseif model != Boussinesq() && iout == 1
            @views file["p"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
            ] = pstrattfc[i0:i1, j0:j1, k0:k1] .* rhoref .* thetaref
        end

        # Write the density fluctuations.
        if prepare_restart || :rhop in output_variables
            HDF5.set_extent_dims(file["rhop"], (sizex, sizey, sizez, iout))
            if model == Boussinesq()
                @views file["rhop"][
                    (io + 1):(io + nx),
                    (jo + 1):(jo + ny),
                    (ko + 1):(ko + nz),
                    iout,
                ] = rhop[i0:i1, j0:j1, k0:k1] .* rhoref
            else
                @views file["rhop"][
                    (io + 1):(io + nx),
                    (jo + 1):(jo + ny),
                    (ko + 1):(ko + nz),
                    iout,
                ] = rho[i0:i1, j0:j1, k0:k1] .* rhoref
            end
        end

        # Write the zonal winds.
        if :u in output_variables
            HDF5.set_extent_dims(file["u"], (sizex, sizey, sizez, iout))
            @views file["u"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] =
                (
                    u[i0:i1, j0:j1, k0:k1] .+
                    u[(i0 - 1):(i1 - 1), j0:j1, k0:k1]
                ) ./ 2 .* uref
        end

        # Write the staggered zonal winds.
        if prepare_restart || :us in output_variables
            HDF5.set_extent_dims(file["us"], (sizex, sizey, sizez, iout))
            @views file["us"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] = u[i0:i1, j0:j1, k0:k1] .* uref
        end

        # Write the meridional winds.
        if :v in output_variables
            HDF5.set_extent_dims(file["v"], (sizex, sizey, sizez, iout))
            @views file["v"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] =
                (
                    v[i0:i1, j0:j1, k0:k1] .+
                    v[i0:i1, (j0 - 1):(j1 - 1), k0:k1]
                ) ./ 2 .* uref
        end

        # Write the staggered meridional winds.
        if prepare_restart || :vs in output_variables
            HDF5.set_extent_dims(file["vs"], (sizex, sizey, sizez, iout))
            @views file["vs"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] = v[i0:i1, j0:j1, k0:k1] .* uref
        end

        # Write the vertical winds.
        if :w in output_variables
            HDF5.set_extent_dims(file["w"], (sizex, sizey, sizez, iout))
            for k in 1:nz, j in 1:ny, i in 1:nx
                file["w"][io + i, jo + j, ko + k, iout] =
                    (
                        compute_vertical_wind(
                            i + i0 - 1,
                            j + j0 - 1,
                            k + k0 - 1,
                            predictands,
                            grid,
                        ) + compute_vertical_wind(
                            i + i0 - 1,
                            j + j0 - 1,
                            k + k0 - 2,
                            predictands,
                            grid,
                        )
                    ) / 2 * uref
            end
        end

        # Write the staggered vertical winds.
        if :ws in output_variables
            HDF5.set_extent_dims(file["ws"], (sizex, sizey, sizez, iout))
            for k in 1:nz, j in 1:ny, i in 1:nx
                file["ws"][io + i, jo + j, ko + k, iout] =
                    compute_vertical_wind(
                        i + i0 - 1,
                        j + j0 - 1,
                        k + k0 - 1,
                        predictands,
                        grid,
                    ) * uref
            end
        end

        # Write the transformed vertical winds.
        if :wtfc in output_variables
            HDF5.set_extent_dims(file["wtfc"], (sizex, sizey, sizez, iout))
            @views file["wtfc"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] =
                (
                    w[i0:i1, j0:j1, k0:k1] .+
                    w[i0:i1, j0:j1, (k0 - 1):(k1 - 1)]
                ) ./ 2 .* uref
        end

        # Write the staggered transformed vertical winds.
        if prepare_restart || :wstfc in output_variables
            HDF5.set_extent_dims(file["wstfc"], (sizex, sizey, sizez, iout))
            @views file["wstfc"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] = w[i0:i1, j0:j1, k0:k1] .* uref
        end

        # Write the potential-temperature fluctuations.
        if :thetap in output_variables
            HDF5.set_extent_dims(file["thetap"], (sizex, sizey, sizez, iout))
            if model == Boussinesq()
                @views file["thetap"][
                    (io + 1):(io + nx),
                    (jo + 1):(jo + ny),
                    (ko + 1):(ko + nz),
                    iout,
                ] =
                    (
                        pstrattfc[i0:i1, j0:j1, k0:k1] ./ (
                            rhostrattfc[i0:i1, j0:j1, k0:k1] .+
                            rhop[i0:i1, j0:j1, k0:k1]
                        ) .- thetastrattfc[i0:i1, j0:j1, k0:k1]
                    ) .* thetaref
            else
                @views file["thetap"][
                    (io + 1):(io + nx),
                    (jo + 1):(jo + ny),
                    (ko + 1):(ko + nz),
                    iout,
                ] =
                    (
                        pstrattfc[i0:i1, j0:j1, k0:k1] ./ (
                            rhostrattfc[i0:i1, j0:j1, k0:k1] .+
                            rho[i0:i1, j0:j1, k0:k1]
                        ) .- thetastrattfc[i0:i1, j0:j1, k0:k1]
                    ) .* thetaref
            end
        end

        # Write the Exner-pressure fluctuations.
        if prepare_restart || :pip in output_variables
            HDF5.set_extent_dims(file["pip"], (sizex, sizey, sizez, iout))
            @views file["pip"][
                (io + 1):(io + nx),
                (jo + 1):(jo + ny),
                (ko + 1):(ko + nz),
                iout,
            ] = pip[i0:i1, j0:j1, k0:k1]
        end

        # Write WKB variables.
        if typeof(testcase) <: AbstractWKBTestCase

            # Write ray-volume properties.
            if prepare_restart || save_ray_volumes
                dk0 = ko == 0 ? 1 : 0
                dk1 = ko + nzz == sizezz ? 1 : 0

                for (output_name, field_name) in zip(
                    ("xr", "yr", "zr", "dxr", "dyr", "dzr"),
                    (:x, :y, :z, :dxray, :dyray, :dzray),
                )
                    HDF5.set_extent_dims(
                        file[output_name],
                        (nray_max, sizex, sizey, sizez + 2, iout),
                    )
                    @views file[output_name][
                        1:nray_max,
                        (io + 1):(io + nx),
                        (jo + 1):(jo + ny),
                        (ko + 2 - dk0):(ko + nz + 1 + dk1),
                        iout,
                    ] =
                        getfield(rays, field_name)[
                            1:nray_max,
                            i0:i1,
                            j0:j1,
                            (k0 - dk0):(k1 + dk1),
                        ] .* lref
                end

                for (output_name, field_name) in zip(
                    ("kr", "lr", "mr", "dkr", "dlr", "dmr"),
                    (:k, :l, :m, :dkray, :dlray, :dmray),
                )
                    HDF5.set_extent_dims(
                        file[output_name],
                        (nray_max, sizex, sizey, sizez + 2, iout),
                    )
                    @views file[output_name][
                        1:nray_max,
                        (io + 1):(io + nx),
                        (jo + 1):(jo + ny),
                        (ko + 2 - dk0):(ko + nz + 1 + dk1),
                        iout,
                    ] =
                        getfield(rays, field_name)[
                            1:nray_max,
                            i0:i1,
                            j0:j1,
                            (k0 - dk0):(k1 + dk1),
                        ] ./ lref
                end

                HDF5.set_extent_dims(
                    file["nr"],
                    (nray_max, sizex, sizey, sizez + 2, iout),
                )
                @views file["nr"][
                    1:nray_max,
                    (io + 1):(io + nx),
                    (jo + 1):(jo + ny),
                    (ko + 2 - dk0):(ko + nz + 1 + dk1),
                    iout,
                ] =
                    rays.dens[
                        1:nray_max,
                        i0:i1,
                        j0:j1,
                        (k0 - dk0):(k1 + dk1),
                    ] .* rhoref .* uref .^ 2 .* tref .* lref .^ dim
            end

            # Write GW tendencies.
            for (field, scaling) in zip(
                (:dudt, :dvdt, :dthetadt),
                (
                    rhoref * uref / tref,
                    rhoref * uref / tref,
                    rhoref * thetaref / tref,
                ),
            )
                if field in output_variables
                    HDF5.set_extent_dims(
                        file[string(field)],
                        (sizex, sizey, sizez, iout),
                    )
                    @views file[string(field)][
                        (io + 1):(io + nx),
                        (jo + 1):(jo + ny),
                        (ko + 1):(ko + nz),
                        iout,
                    ] =
                        getfield(tendencies, field)[i0:i1, j0:j1, k0:k1] .*
                        scaling
                end
            end
        end

        # Return.
        return
    end

    # Return.
    return iout
end
