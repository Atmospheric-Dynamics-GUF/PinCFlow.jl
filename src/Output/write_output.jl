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

The output is written in parallel, using the chunking prepared by `create_output`. The grid, i.e. the fields `x`, `y` and `zc` of `state.grid`, as well as the fields of `state.atmosphere` are only written if `iout == 1` (which should only be the case for the initial output). In Boussinesq mode, the fields of `state.atmosphere` do not have a spatial dependence and are therefore not written at all. In compressible mode, the mass-weighted potential temperature and squared buoyancy frequency have a temporal dependence and are therefore written even if `iout != 1`. Any other field is only written if it is listed in `state.namelists.output.output_variables` or if it is essential for restarts and `state.namelists.output.prepare_restart == true`.

The list of available output variables (as specified in `state.namelists.output.output_variables`) is as follows.

  - `:rhop`: Density fluctuations (restart variable).

  - `:u`: Zonal wind.

  - `:us`: Staggered zonal wind (restart variable).

  - `:v`: Meridional wind.

  - `:vs`: Staggered meridional wind (restart variable).

  - `:w`: Vertical wind (computed with `compute_vertical_wind`).

  - `:ws`: Staggered vertical wind (computed with `compute_vertical_wind`).

  - `:wt`: Transformed vertical wind.

  - `:wts`: Staggered transformed vertical wind (restart variable).

  - `:thetap`: Potential-temperature fluctuations.

  - `:pip`: Exner-pressure fluctuations (restart variable).

  - `:dudt`: Zonal-momentum drag due to unresolved gravity waves.

  - `:dvdt`: Meridional-momentum drag due to unresolved gravity waves.

  - `:dthetadt`: Mass-weighted potential-temperature tendency due to unresolved gravity waves.

  - `:dchidt`: Leading-order tracer impact of unresolved gravity waves.

  - `:uchi`: Zonal tracer fluxes due to unresolved gravity waves.

  - `:vchi`: Meridional tracer fluxes due to unresolved gravity waves.

  - `:wchi`: Vertical tracer fluxes due to unresolved gravity waves.

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
    (; domain, grid) = state
    (; x_size, y_size, z_size) = state.namelists.domain
    (; prepare_restart, save_ray_volumes, output_variables, output_file) =
        state.namelists.output
    (; model) = state.namelists.atmosphere
    (; wkb_mode) = state.namelists.wkb
    (; comm, master, nx, ny, nz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; tref, lref, rhoref, thetaref, uref) = state.constants
    (; x, y, zc, zctilde) = grid
    (; rhobar, thetabar, n2, pbar) = state.atmosphere
    (; predictands) = state.variables
    (; rho, rhop, u, v, w, pip, p) = predictands
    (; nray_max, rays, tendencies) = state.wkb

    # Print information.
    if master
        println(repeat("-", 80))
        println("Output into file ", output_file)
        println("Physical time: ", time * tref, " s")
        println("Machine time: ", canonicalize(now() - machine_start_time))
        println(repeat("-", 80))
        println("")
    end

    # Advance output counter.
    iout += 1

    # Determine dimensionality.
    dim = 1
    if x_size > 1
        dim += 1
    end
    if y_size > 1
        dim += 1
    end

    # Define slices.
    dk0 = ko == 0 ? 1 : 0
    (rr, ii, jj, kk, kkr) = (1:nray_max, i0:i1, j0:j1, k0:k1, (k0 - dk0):k1)
    (iid, jjd, kkd, kkrd) = (
        (io + 1):(io + nx),
        (jo + 1):(jo + ny),
        (ko + 1):(ko + nz),
        (ko + 2 - dk0):(ko + nz + 1),
    )

    # Open the file. Note: Fused in-place assignments cannot be used here!
    @ivy h5open(output_file, "r+", comm) do file

        # Write the time.
        HDF5.set_extent_dims(file["t"], (iout,))
        file["t"][iout] = time * tref

        # Write the horizontal grid.
        if iout == 1
            file["x"][iid] = x[ii] .* lref
            file["y"][jjd] = y[jj] .* lref
        end

        # Write the vertical grid.
        if iout == 1
            file["z"][iid, jjd, kkd] = zc[ii, jj, kk] .* lref
            file["ztilde"][iid, jjd, kkrd] = zctilde[ii, jj, kkr] .* lref
        end

        # Write the background density.
        if model != Boussinesq() && iout == 1
            file["rhobar"][iid, jjd, kkd] = rhobar[ii, jj, kk] .* rhoref
        end

        # Write the background potential temperature.
        if model != Boussinesq() && iout == 1
            file["thetabar"][iid, jjd, kkd] = thetabar[ii, jj, kk] .* thetaref
        end

        # Write the squared buoyancy frequency.
        if model != Boussinesq() && iout == 1
            file["n2"][iid, jjd, kkd] = n2[ii, jj, kk] ./ tref .^ 2
        end

        # Write the mass-weighted potential temperature.
        if model == Compressible()
            HDF5.set_extent_dims(file["p"], (x_size, y_size, z_size, iout))
            file["p"][iid, jjd, kkd, iout] = p[ii, jj, kk] .* rhoref .* thetaref
        elseif model != Boussinesq() && iout == 1
            file["p"][iid, jjd, kkd] = pbar[ii, jj, kk] .* rhoref .* thetaref
        end

        # Write the density fluctuations.
        if prepare_restart || :rhop in output_variables
            HDF5.set_extent_dims(file["rhop"], (x_size, y_size, z_size, iout))
            if model == Boussinesq()
                file["rhop"][iid, jjd, kkd, iout] = rhop[ii, jj, kk] .* rhoref
            else
                file["rhop"][iid, jjd, kkd, iout] = rho[ii, jj, kk] .* rhoref
            end
        end

        # Write the zonal winds.
        if :u in output_variables
            HDF5.set_extent_dims(file["u"], (x_size, y_size, z_size, iout))
            file["u"][iid, jjd, kkd, iout] =
                map(CartesianIndices((ii, jj, kk))) do ijk
                    (i, j, k) = Tuple(ijk)
                    return (u[i, j, k] + u[i - 1, j, k]) / 2 * uref
                end
        end

        # Write the staggered zonal winds.
        if prepare_restart || :us in output_variables
            HDF5.set_extent_dims(file["us"], (x_size, y_size, z_size, iout))
            file["us"][iid, jjd, kkd, iout] = u[ii, jj, kk] .* uref
        end

        # Write the meridional winds.
        if :v in output_variables
            HDF5.set_extent_dims(file["v"], (x_size, y_size, z_size, iout))
            file["v"][iid, jjd, kkd, iout] =
                map(CartesianIndices((ii, jj, kk))) do ijk
                    (i, j, k) = Tuple(ijk)
                    return (v[i, j, k] + v[i, j - 1, k]) / 2 * uref
                end
        end

        # Write the staggered meridional winds.
        if prepare_restart || :vs in output_variables
            HDF5.set_extent_dims(file["vs"], (x_size, y_size, z_size, iout))
            file["vs"][iid, jjd, kkd, iout] = v[ii, jj, kk] .* uref
        end

        # Write the vertical winds.
        if :w in output_variables
            HDF5.set_extent_dims(file["w"], (x_size, y_size, z_size, iout))
            file["w"][iid, jjd, kkd, iout] =
                map(CartesianIndices((ii, jj, kk))) do ijk
                    (i, j, k) = Tuple(ijk)
                    return (
                        compute_vertical_wind(i, j, k, state) +
                        compute_vertical_wind(i, j, k - 1, state)
                    ) / 2 * uref
                end
        end

        # Write the staggered vertical winds.
        if :ws in output_variables
            HDF5.set_extent_dims(file["ws"], (x_size, y_size, z_size, iout))
            file["ws"][iid, jjd, kkd, iout] =
                map(CartesianIndices((ii, jj, kk))) do ijk
                    (i, j, k) = Tuple(ijk)
                    return compute_vertical_wind(i, j, k, state) * uref
                end
        end

        # Write the transformed vertical winds.
        if :wt in output_variables
            HDF5.set_extent_dims(file["wt"], (x_size, y_size, z_size, iout))
            file["wt"][iid, jjd, kkd, iout] =
                map(CartesianIndices((ii, jj, kk))) do ijk
                    (i, j, k) = Tuple(ijk)
                    return (w[i, j, k] + w[i, j, k - 1]) / 2 * uref
                end
        end

        # Write the staggered transformed vertical winds.
        if prepare_restart || :wts in output_variables
            HDF5.set_extent_dims(file["wts"], (x_size, y_size, z_size, iout))
            file["wts"][iid, jjd, kkd, iout] = w[ii, jj, kk] .* uref
        end

        # Write the potential-temperature fluctuations.
        if :thetap in output_variables
            HDF5.set_extent_dims(file["thetap"], (x_size, y_size, z_size, iout))
            if model == Boussinesq()
                file["thetap"][iid, jjd, kkd, iout] =
                    (
                        pbar[ii, jj, kk] ./
                        (rhobar[ii, jj, kk] .+ rhop[ii, jj, kk]) .-
                        thetabar[ii, jj, kk]
                    ) .* thetaref
            else
                file["thetap"][iid, jjd, kkd, iout] =
                    (
                        pbar[ii, jj, kk] ./
                        (rhobar[ii, jj, kk] .+ rho[ii, jj, kk]) .-
                        thetabar[ii, jj, kk]
                    ) .* thetaref
            end
        end

        # Write the Exner-pressure fluctuations.
        if prepare_restart || :pip in output_variables
            HDF5.set_extent_dims(file["pip"], (x_size, y_size, z_size, iout))
            file["pip"][iid, jjd, kkd, iout] = pip[ii, jj, kk]
        end

        if !(typeof(state.namelists.tracer.tracer_setup) <: NoTracer)
            for field in fieldnames(TracerPredictands)
                HDF5.set_extent_dims(
                    file[string(field)],
                    (x_size, y_size, z_size, iout),
                )
                file[string(field)][iid, jjd, kkd, iout] =
                    getfield(state.tracer.tracerpredictands, field)[
                        ii,
                        jj,
                        kk,
                    ] ./ (rhobar[ii, jj, kk] .+ rho[ii, jj, kk])
            end

            if state.namelists.tracer.leading_order_impact &&
               :dchidt in output_variables
                for field in (:dchidt,)
                    HDF5.set_extent_dims(
                        file[string(field)],
                        (x_size, y_size, z_size, iout),
                    )
                    @views file[string(field)][iid, jjd, kkd, iout] =
                        getfield(state.tracer.tracerforcings.chiq0, field)[
                            ii,
                            jj,
                            kk,
                        ] ./ tref ./ (rhobar[ii, jj, kk] .+ rho[ii, jj, kk])
                end
                for field in (:uchi, :vchi, :wchi)
                    HDF5.set_extent_dims(
                        file[string(field)],
                        (x_size, y_size, z_size, iout),
                    )
                    @views file[string(field)][iid, jjd, kkd, iout] =
                        getfield(state.tracer.tracerforcings.chiq0, field)[
                            ii,
                            jj,
                            kk,
                        ] .* uref ./ rhobar[ii, jj, kk]
                end
            end
        end

        # Write WKB variables.
        if wkb_mode != NoWKB()

            # Write ray-volume properties.
            if prepare_restart || save_ray_volumes
                for (output_name, field_name) in zip(
                    ("xr", "yr", "zr", "dxr", "dyr", "dzr"),
                    (:x, :y, :z, :dxray, :dyray, :dzray),
                )
                    HDF5.set_extent_dims(
                        file[output_name],
                        (nray_max, x_size, y_size, z_size + 1, iout),
                    )
                    file[output_name][1:nray_max, iid, jjd, kkrd, iout] =
                        getfield(rays, field_name)[rr, ii, jj, kkr] .* lref
                end

                for (output_name, field_name) in zip(
                    ("kr", "lr", "mr", "dkr", "dlr", "dmr"),
                    (:k, :l, :m, :dkray, :dlray, :dmray),
                )
                    HDF5.set_extent_dims(
                        file[output_name],
                        (nray_max, x_size, y_size, z_size + 1, iout),
                    )
                    file[output_name][1:nray_max, iid, jjd, kkrd, iout] =
                        getfield(rays, field_name)[rr, ii, jj, kkr] ./ lref
                end

                HDF5.set_extent_dims(
                    file["nr"],
                    (nray_max, x_size, y_size, z_size + 1, iout),
                )
                file["nr"][1:nray_max, iid, jjd, kkrd, iout] =
                    rays.dens[rr, ii, jj, kkr] .* rhoref .* uref .^ 2 .* tref .*
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
                        (x_size, y_size, z_size, iout),
                    )
                    file[string(field)][iid, jjd, kkd, iout] =
                        getfield(tendencies, field)[ii, jj, kk] .* scaling
                end
            end
        end

        return
    end

    return iout
end
