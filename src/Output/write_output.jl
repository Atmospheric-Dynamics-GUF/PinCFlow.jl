"""
```julia
write_output(
    state::State,
    time::AbstractFloat,
    iout::Integer,
    machine_start_time::DateTime,
)::Integer
```

Write the current simulation state to a previously created HDF5 output file and return the advanced output counter `iout`.

The output is written in parallel, using the chunking prepared by `create_output`. The grid, i.e. the fields `x`, `y` and `ztfc` of `state.grid`, as well as the fields of `state.atmosphere` are only written if `iout == 1` (which should only be the case for the initial output). In Boussinesq mode, the fields of `state.atmosphere` do not have a spatial dependence and are therefore not written at all. In compressible mode, the mass-weighted potential temperature and squared buoyancy frequency have a temporal dependence and are therefore written even if `iout != 1`. Any other field is only written if it is listed in `state.namelists.output.output_variables` or if it is essential for restarts and `state.namelists.output.prepare_restart == true`.

The list of available output variables (as specified in `state.namelists.output.output_variables`) is as follows.

  - `:rhop`: Density fluctuations (restart variable).

  - `:u`: Zonal wind.

  - `:us`: Staggered zonal wind (restart variable).

  - `:v`: Meridional wind.

  - `:vs`: Staggered meridional wind (restart variable).

  - `:w`: Vertical wind (computed with `compute_vertical_wind`).

  - `:ws`: Staggered vertical wind (computed with `compute_vertical_wind`).

  - `:wtfc`: Transformed vertical wind.

  - `:wstfc`: Staggered transformed vertical wind (restart variable).

  - `:thetap`: Potential-temperature fluctuations.

  - `:pip`: Exner-pressure fluctuations (restart variable).

  - `:dudt`: Zonal-momentum drag due to unresolved gravity waves.

  - `:dvdt`: Meridional-momentum drag due to unresolved gravity waves.

  - `:dthetadt`: Mass-weighted potential-temperature tendency due to unresolved gravity waves.

An output of all ray-volume properties is provided if `state.namelists.output.save_ray_volumes == true` and/or `state.namelists.output.prepare_restart == true`.

All output variables are re-dimensionalized with the scale parameters stored in `state.constants`.

# Arguments

  - `state`: Model state.

  - `time`: Simulation time.

  - `iout`: Output counter. This is the temporal index of the output. It is advanced before the output is written, so that the first call of `write_output` should receive `iout = 0`.

  - `machine_start_time`: Wall-clock start time.

# See also

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)
"""
function write_output end

function write_output(
    state::State,
    time::AbstractFloat,
    iout::Integer,
    machine_start_time::DateTime,
)::Integer

    # Get all necessary fields.
    (; domain, grid) = state
    (; sizex, sizey, sizez) = state.namelists.domain
    (; prepare_restart, save_ray_volumes, output_variables, output_file) =
        state.namelists.output
    (; model, testcase) = state.namelists.setting
    (; comm, master, nx, ny, nz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
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

    # Define slices.
    dk0 = ko == 0 ? 1 : 0
    (r, i, j, k, kr) = (1:nray_max, i0:i1, j0:j1, k0:k1, (k0 - dk0):k1)
    (im1, jm1, km1) = ((i0 - 1):(i1 - 1), (j0 - 1):(j1 - 1), (k0 - 1):(k1 - 1))
    (id, jd, kd, krd) = (
        (io + 1):(io + nx),
        (jo + 1):(jo + ny),
        (ko + 1):(ko + nz),
        (ko + 2 - dk0):(ko + nz + 1),
    )

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
            @views file["z"][id, jd, kd] = ztfc[i, j, k] .* lref
        end

        # Write the background density.
        if model != Boussinesq() && iout == 1
            @views file["rhobar"][id, jd, kd] = rhostrattfc[i, j, k] .* rhoref
        end

        # Write the background potential temperature.
        if model != Boussinesq() && iout == 1
            @views file["thetabar"][id, jd, kd] =
                thetastrattfc[i, j, k] .* thetaref
        end

        # Write the squared buoyancy frequency.
        if model != Boussinesq() && iout == 1
            @views file["n2"][id, jd, kd] = bvsstrattfc[i, j, k] ./ tref .^ 2
        end

        # Write the mass-weighted potential temperature.
        if model == Compressible()
            HDF5.set_extent_dims(file["p"], (sizex, sizey, sizez, iout))
            @views file["p"][id, jd, kd, iout] =
                p[i, j, k] .* rhoref .* thetaref
        elseif model != Boussinesq() && iout == 1
            @views file["p"][id, jd, kd] =
                pstrattfc[i, j, k] .* rhoref .* thetaref
        end

        # Write the density fluctuations.
        if prepare_restart || :rhop in output_variables
            HDF5.set_extent_dims(file["rhop"], (sizex, sizey, sizez, iout))
            if model == Boussinesq()
                @views file["rhop"][id, jd, kd, iout] = rhop[i, j, k] .* rhoref
            else
                @views file["rhop"][id, jd, kd, iout] = rho[i, j, k] .* rhoref
            end
        end

        # Write the zonal winds.
        if :u in output_variables
            HDF5.set_extent_dims(file["u"], (sizex, sizey, sizez, iout))
            @views file["u"][id, jd, kd, iout] =
                (u[i, j, k] .+ u[im1, j, k]) ./ 2 .* uref
        end

        # Write the staggered zonal winds.
        if prepare_restart || :us in output_variables
            HDF5.set_extent_dims(file["us"], (sizex, sizey, sizez, iout))
            @views file["us"][id, jd, kd, iout] = u[i, j, k] .* uref
        end

        # Write the meridional winds.
        if :v in output_variables
            HDF5.set_extent_dims(file["v"], (sizex, sizey, sizez, iout))
            @views file["v"][id, jd, kd, iout] =
                (v[i, j, k] .+ v[i, jm1, k]) ./ 2 .* uref
        end

        # Write the staggered meridional winds.
        if prepare_restart || :vs in output_variables
            HDF5.set_extent_dims(file["vs"], (sizex, sizey, sizez, iout))
            @views file["vs"][id, jd, kd, iout] = v[i, j, k] .* uref
        end

        # Write the vertical winds.
        if :w in output_variables
            HDF5.set_extent_dims(file["w"], (sizex, sizey, sizez, iout))
            file["w"][id, jd, kd, iout] = (
                ijk ->
                    (
                        compute_vertical_wind(
                            ijk[1],
                            ijk[2],
                            ijk[3],
                            predictands,
                            grid,
                        ) + compute_vertical_wind(
                            ijk[1],
                            ijk[2],
                            ijk[3] - 1,
                            predictands,
                            grid,
                        )
                    ) / 2 * uref
            ).(
                CartesianIndices((i, j, k)),
            )
        end

        # Write the staggered vertical winds.
        if :ws in output_variables
            HDF5.set_extent_dims(file["ws"], (sizex, sizey, sizez, iout))
            file["ws"][id, jd, kd, iout] = (
                ijk ->
                    compute_vertical_wind(
                        ijk[1],
                        ijk[2],
                        ijk[3],
                        predictands,
                        grid,
                    ) * uref
            ).(
                CartesianIndices((i, j, k)),
            )
        end

        # Write the transformed vertical winds.
        if :wtfc in output_variables
            HDF5.set_extent_dims(file["wtfc"], (sizex, sizey, sizez, iout))
            @views file["wtfc"][id, jd, kd, iout] =
                (w[i, j, k] .+ w[i, j, km1]) ./ 2 .* uref
        end

        # Write the staggered transformed vertical winds.
        if prepare_restart || :wstfc in output_variables
            HDF5.set_extent_dims(file["wstfc"], (sizex, sizey, sizez, iout))
            @views file["wstfc"][id, jd, kd, iout] = w[i, j, k] .* uref
        end

        # Write the potential-temperature fluctuations.
        if :thetap in output_variables
            HDF5.set_extent_dims(file["thetap"], (sizex, sizey, sizez, iout))
            if model == Boussinesq()
                @views file["thetap"][id, jd, kd, iout] =
                    (
                        pstrattfc[i, j, k] ./
                        (rhostrattfc[i, j, k] .+ rhop[i, j, k]) .-
                        thetastrattfc[i, j, k]
                    ) .* thetaref
            else
                @views file["thetap"][id, jd, kd, iout] =
                    (
                        pstrattfc[i, j, k] ./
                        (rhostrattfc[i, j, k] .+ rho[i, j, k]) .-
                        thetastrattfc[i, j, k]
                    ) .* thetaref
            end
        end

        # Write the Exner-pressure fluctuations.
        if prepare_restart || :pip in output_variables
            HDF5.set_extent_dims(file["pip"], (sizex, sizey, sizez, iout))
            @views file["pip"][id, jd, kd, iout] = pip[i, j, k]
        end

        if !(typeof(state.namelists.tracer.tracersetup) <: NoTracer)
            for field in fieldnames(TracerPredictands)
                HDF5.set_extent_dims(
                    file[string(field)],
                    (sizex, sizey, sizez, iout),
                )
                @views file[string(field)][id, jd, kd, iout] =
                    getfield(state.tracer.tracerpredictands, field)[i, j, k] ./
                    (rhostrattfc[i, j, k] .+ rho[i, j, k]) .* lref
            end
        end

        if !(typeof(state.namelists.ice.icesetup) <: NoIce)
            for field in fieldnames(IcePredictands)
                HDF5.set_extent_dims(
                    file[string(field)],
                    (sizex, sizey, sizez, iout),
                )
                @views file[string(field)][id, jd, kd, iout] =
                    getfield(state.ice.icepredictands, field)[i, j, k] ./
                    (rhostrattfc[i, j, k] .+ rho[i, j, k]) .* lref
            end
        end

        if !(typeof(state.namelists.turbulence.turbulencesetup) <: NoTurbulence)
            for field in fieldnames(TurbulencePredictands)
                HDF5.set_extent_dims(
                    file[string(field)],
                    (sizex, sizey, sizez, iout),
                )
                @views file[string(field)][id, jd, kd, iout] =
                    getfield(state.turbulence.turbulencepredictands, field)[
                        i,
                        j,
                        k,
                    ] ./ (rhostrattfc[i, j, k] .+ rho[i, j, k]) .* lref .^ 2 ./
                    tref .^ 2
            end
        end

        # Write WKB variables.
        if typeof(testcase) <: AbstractWKBTestCase

            # Write ray-volume properties.
            if prepare_restart || save_ray_volumes
                for (output_name, field_name) in zip(
                    ("xr", "yr", "zr", "dxr", "dyr", "dzr"),
                    (:x, :y, :z, :dxray, :dyray, :dzray),
                )
                    HDF5.set_extent_dims(
                        file[output_name],
                        (nray_max, sizex, sizey, sizez + 1, iout),
                    )
                    @views file[output_name][1:nray_max, id, jd, krd, iout] =
                        getfield(rays, field_name)[r, i, j, kr] .* lref
                end

                for (output_name, field_name) in zip(
                    ("kr", "lr", "mr", "dkr", "dlr", "dmr"),
                    (:k, :l, :m, :dkray, :dlray, :dmray),
                )
                    HDF5.set_extent_dims(
                        file[output_name],
                        (nray_max, sizex, sizey, sizez + 1, iout),
                    )
                    @views file[output_name][1:nray_max, id, jd, krd, iout] =
                        getfield(rays, field_name)[r, i, j, kr] ./ lref
                end

                HDF5.set_extent_dims(
                    file["nr"],
                    (nray_max, sizex, sizey, sizez + 1, iout),
                )
                @views file["nr"][1:nray_max, id, jd, krd, iout] =
                    rays.dens[r, i, j, kr] .* rhoref .* uref .^ 2 .* tref .*
                    lref .^ dim
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
                    @views file[string(field)][id, jd, kd, iout] =
                        getfield(tendencies, field)[i, j, k] .* scaling
                end
            end
        end

        # Return.
        return
    end

    # Return.
    return iout
end
