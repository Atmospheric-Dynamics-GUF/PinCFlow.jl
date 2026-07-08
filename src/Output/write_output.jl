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

  - `:uu`: Zonal zonal-momentum flux due to unresolved gravity waves.

  - `:uv`: Zonal meridional-momentum flux due to unresolved gravity waves.

  - `:uw`: Zonal vertical-momentum flux due to unresolved gravity waves.

  - `:vv`: Meridional meridional-momentum flux due to unresolved gravity waves.

  - `:vw`: Meridional vertical-momentum flux due to unresolved gravity waves.

  - `:utheta`: Zonal mass-weighted potential-temperature flux due to unresolved gravity waves.

  - `:vtheta`: Meridional mass-weighted potential-temperature flux due to unresolved gravity waves.

  - `:e`: Energy density of unresolved gravity waves.

  - `:dudt`: Zonal-momentum drag due to unresolved gravity waves.

  - `:dvdt`: Meridional-momentum drag due to unresolved gravity waves.

  - `:dthetadt`: Mass-weighted potential-temperature tendency due to unresolved gravity waves.

  - `:dchidt0`: Leading-order tracer impact of unresolved gravity waves.

  - `:uchi0`: Zonal tracer fluxes due to unresolved gravity waves.

  - `:vchi0`: Meridional tracer fluxes due to unresolved gravity waves.

  - `:wchi0`: Vertical tracer fluxes due to unresolved gravity waves.

  - `:tke`: Turbulent kinetic energy.

  - `:shear_production`: Turbulent kinetic energy production due to wind shear.

  - `:buoyancy_production`: Turbulent kinetic energy production/destruction due to buoyancy.

  - `:launch_mode_count`: Numbers of modes selected by the elastic-mode-selection algorithm.

  - `:launch_power_fraction`: Power fractions retained by the elastic-mode-selection algorithm.

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

@ivy function write_output(
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
    (; wkb_mode, elastic_mode_selection) = state.namelists.wkb
    (; comm, master, nx, ny, nz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; tref, lref, rhoref, thetaref, uref) = state.constants
    (; x, y, zc, zctilde) = grid
    (; rhobar, thetabar, n2, pbar) = state.atmosphere
    (; predictands) = state.variables
    (; rho, rhop, u, v, w, pip, p) = predictands
    (; bins, rays, tendencies, integrals) = state.wkb

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
    (rr, ii, jj, kk, kkr) = (1:bins, i0:i1, j0:j1, k0:k1, (k0 - dk0):k1)
    (iid, jjd, kkd, kkrd) = (
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
            file["x"][iid] = @share x[ii] * lref
            file["y"][jjd] = @share y[jj] * lref
        end

        # Write the vertical grid.
        if iout == 1
            file["z"][iid, jjd, kkd] = @share zc[ii, jj, kk] * lref
            file["ztilde"][iid, jjd, kkrd] = @share zctilde[ii, jj, kkr] * lref
        end

        # Write the background density.
        if model !== :Boussinesq && iout == 1
            file["rhobar"][iid, jjd, kkd] = @share rhobar[ii, jj, kk] * rhoref
        end

        # Write the background potential temperature.
        if model !== :Boussinesq && iout == 1
            file["thetabar"][iid, jjd, kkd] =
                @share thetabar[ii, jj, kk] * thetaref
        end

        # Write the squared buoyancy frequency.
        if model !== :Boussinesq && iout == 1
            file["n2"][iid, jjd, kkd] = @share n2[ii, jj, kk] / tref^2
        end

        # Write the mass-weighted potential temperature.
        if model === :Compressible
            HDF5.set_extent_dims(file["p"], (x_size, y_size, z_size, iout))
            file["p"][iid, jjd, kkd, iout] =
                @share p[ii, jj, kk] * rhoref * thetaref
        elseif model !== :Boussinesq && iout == 1
            file["p"][iid, jjd, kkd] =
                @share pbar[ii, jj, kk] * rhoref * thetaref
        end

        # Write the density fluctuations.
        if prepare_restart || :rhop in output_variables
            HDF5.set_extent_dims(file["rhop"], (x_size, y_size, z_size, iout))
            if model === :Boussinesq
                file["rhop"][iid, jjd, kkd, iout] =
                    @share rhop[ii, jj, kk] * rhoref
            else
                file["rhop"][iid, jjd, kkd, iout] =
                    @share rho[ii, jj, kk] * rhoref
            end
        end

        # Write the zonal winds.
        if :u in output_variables
            uc = zeros(nx, ny, nz)
            @share for k in kk, j in jj, i in ii
                uc[i - i0 + 1, j - j0 + 1, k - k0 + 1] =
                    (u[i, j, k] + u[i - 1, j, k]) / 2 * uref
            end

            HDF5.set_extent_dims(file["u"], (x_size, y_size, z_size, iout))
            file["u"][iid, jjd, kkd, iout] = uc
        end

        # Write the staggered zonal winds.
        if prepare_restart || :us in output_variables
            HDF5.set_extent_dims(file["us"], (x_size, y_size, z_size, iout))
            file["us"][iid, jjd, kkd, iout] = @share u[ii, jj, kk] * uref
        end

        # Write the meridional winds.
        if :v in output_variables
            vc = zeros(nx, ny, nz)
            @share for k in kk, j in jj, i in ii
                vc[i - i0 + 1, j - j0 + 1, k - k0 + 1] =
                    (v[i, j, k] + v[i, j - 1, k]) / 2 * uref
            end

            HDF5.set_extent_dims(file["v"], (x_size, y_size, z_size, iout))
            file["v"][iid, jjd, kkd, iout] = vc
        end

        # Write the staggered meridional winds.
        if prepare_restart || :vs in output_variables
            HDF5.set_extent_dims(file["vs"], (x_size, y_size, z_size, iout))
            file["vs"][iid, jjd, kkd, iout] = @share v[ii, jj, kk] * uref
        end

        # Write the vertical winds.
        if :w in output_variables
            wc = zeros(nx, ny, nz)
            @share for k in kk, j in jj, i in ii
                wc[i - i0 + 1, j - j0 + 1, k - k0 + 1] =
                    (
                        compute_vertical_wind(i, j, k, state) +
                        compute_vertical_wind(i, j, k - 1, state)
                    ) / 2 * uref
            end

            HDF5.set_extent_dims(file["w"], (x_size, y_size, z_size, iout))
            file["w"][iid, jjd, kkd, iout] = wc
        end

        # Write the staggered vertical winds.
        if :ws in output_variables
            wu = zeros(nx, ny, nz)
            @share for k in kk, j in jj, i in ii
                wu[i - i0 + 1, j - j0 + 1, k - k0 + 1] =
                    compute_vertical_wind(i, j, k, state) * uref
            end

            HDF5.set_extent_dims(file["ws"], (x_size, y_size, z_size, iout))
            file["ws"][iid, jjd, kkd, iout] = wu
        end

        # Write the transformed vertical winds.
        if :wt in output_variables
            wtc = zeros(nx, ny, nz)
            @share for k in kk, j in jj, i in ii
                wtc[i - i0 + 1, j - j0 + 1, k - k0 + 1] =
                    (w[i, j, k] + w[i, j, k - 1]) / 2 * uref
            end

            HDF5.set_extent_dims(file["wt"], (x_size, y_size, z_size, iout))
            file["wt"][iid, jjd, kkd, iout] = wtc
        end

        # Write the staggered transformed vertical winds.
        if prepare_restart || :wts in output_variables
            HDF5.set_extent_dims(file["wts"], (x_size, y_size, z_size, iout))
            file["wts"][iid, jjd, kkd, iout] = @share w[ii, jj, kk] * uref
        end

        # Write the potential-temperature fluctuations.
        if :thetap in output_variables
            HDF5.set_extent_dims(file["thetap"], (x_size, y_size, z_size, iout))
            if model === :Boussinesq
                file["thetap"][iid, jjd, kkd, iout] = @share (
                    pbar[ii, jj, kk] / (rhobar[ii, jj, kk] + rhop[ii, jj, kk]) -
                    thetabar[ii, jj, kk]
                ) * thetaref
            else
                file["thetap"][iid, jjd, kkd, iout] = @share (
                    pbar[ii, jj, kk] / (rhobar[ii, jj, kk] + rho[ii, jj, kk]) -
                    thetabar[ii, jj, kk]
                ) * thetaref
            end
        end

        # Write the Exner-pressure fluctuations.
        if prepare_restart || :pip in output_variables
            HDF5.set_extent_dims(file["pip"], (x_size, y_size, z_size, iout))
            file["pip"][iid, jjd, kkd, iout] = pip[ii, jj, kk]
        end

        if state.namelists.tracer.tracer_setup !== :NoTracer
            for field_name in fieldnames(TracerPredictands)
                field = getfield(state.tracer.tracerpredictands, field_name)

                HDF5.set_extent_dims(
                    file[string(field_name)],
                    (x_size, y_size, z_size, iout),
                )
                file[string(field_name)][iid, jjd, kkd, iout] =
                    @share field[ii, jj, kk] /
                           (rhobar[ii, jj, kk] + rho[ii, jj, kk])
            end

            if state.namelists.tracer.leading_order_impact &&
               :dchidt0 in output_variables
                HDF5.set_extent_dims(
                    file["dchidt0"],
                    (x_size, y_size, z_size, iout),
                )
                file["dchidt0"][iid, jjd, kkd, iout] =
                    @share state.tracer.tracerwkbtendencies.dchidt0[
                        ii,
                        jj,
                        kk,
                    ] / tref / (rhobar[ii, jj, kk] + rho[ii, jj, kk])
            end

            if state.namelists.tracer.leading_order_impact &&
               :uchi0 in output_variables
                HDF5.set_extent_dims(
                    file["uchi0"],
                    (x_size, y_size, z_size, iout),
                )
                file["uchi0"][iid, jjd, kkd, iout] =
                    @share state.tracer.tracerwkbintegrals.uchi0[ii, jj, kk] *
                           uref / rhobar[ii, jj, kk]
            end

            if state.namelists.tracer.leading_order_impact &&
               :vchi0 in output_variables
                HDF5.set_extent_dims(
                    file["vchi0"],
                    (x_size, y_size, z_size, iout),
                )
                file["vchi0"][iid, jjd, kkd, iout] =
                    @share state.tracer.tracerwkbintegrals.vchi0[ii, jj, kk] *
                           uref / rhobar[ii, jj, kk]
            end

            if state.namelists.tracer.leading_order_impact &&
               :wchi0 in output_variables
                HDF5.set_extent_dims(
                    file["wchi0"],
                    (x_size, y_size, z_size, iout),
                )
                file["wchi0"][iid, jjd, kkd, iout] =
                    @share state.tracer.tracerwkbintegrals.wchi0[ii, jj, kk] *
                           uref / rhobar[ii, jj, kk]
            end
        end

        if state.namelists.turbulence.turbulence_scheme !== :NoTurbulence
            if prepare_restart || :tke in output_variables
                HDF5.set_extent_dims(
                    file["tke"],
                    (x_size, y_size, z_size, iout),
                )
                file["tke"][iid, jjd, kkd, iout] =
                    @share state.turbulence.turbulencepredictands.tke[
                        ii,
                        jj,
                        kk,
                    ] / (rhobar[ii, jj, kk] + rho[ii, jj, kk]) * (lref^2.0) /
                           (tref^2.0)
            end

            if :shear_production in output_variables
                HDF5.set_extent_dims(
                    file["shear_production"],
                    (x_size, y_size, z_size, iout),
                )
                file["shear_production"][iid, jjd, kkd, iout] =
                    @share state.turbulence.turbulenceauxiliaries.shear_production[
                        ii,
                        jj,
                        kk,
                    ] * uref^2 / tref
            end

            if :buoyancy_production in output_variables
                HDF5.set_extent_dims(
                    file["buoyancy_production"],
                    (x_size, y_size, z_size, iout),
                )
                file["buoyancy_production"][iid, jjd, kkd, iout] =
                    @share state.turbulence.turbulenceauxiliaries.buoyancy_production[
                        ii,
                        jj,
                        kk,
                    ] * uref^2 / tref
            end
        end

        # Write WKB variables.
        if wkb_mode !== :NoWKB

            # Write ray-volume properties.
            if prepare_restart || save_ray_volumes
                for (output_name, field_name) in zip(
                    ("xr", "yr", "zr", "dxr", "dyr", "dzr"),
                    (:x, :y, :z, :dxray, :dyray, :dzray),
                )
                    field = getproperty(rays, field_name)

                    HDF5.set_extent_dims(
                        file[output_name],
                        (bins, x_size, y_size, z_size + 1, iout),
                    )
                    file[output_name][1:bins, iid, jjd, kkrd, iout] =
                        @share field[rr, ii, jj, kkr] * lref
                end

                for (output_name, field_name) in zip(
                    ("kr", "lr", "mr", "dkr", "dlr", "dmr"),
                    (:k, :l, :m, :dkray, :dlray, :dmray),
                )
                    field = getproperty(rays, field_name)

                    HDF5.set_extent_dims(
                        file[output_name],
                        (bins, x_size, y_size, z_size + 1, iout),
                    )
                    file[output_name][1:bins, iid, jjd, kkrd, iout] =
                        @share field[rr, ii, jj, kkr] / lref
                end

                HDF5.set_extent_dims(
                    file["nr"],
                    (bins, x_size, y_size, z_size + 1, iout),
                )
                file["nr"][1:bins, iid, jjd, kkrd, iout] =
                    @share rays.dens[rr, ii, jj, kkr] *
                           rhoref *
                           uref^2 *
                           tref *
                           lref^dim
            end

            # Write GW integrals.
            for (field_name, scaling) in zip(
                (:uu, :uv, :uw, :vv, :vw, :utheta, :vtheta, :e),
                (
                    (rhoref * uref^2 for index in 1:5)...,
                    rhoref * uref * thetaref,
                    rhoref * uref * thetaref,
                    rhoref * uref^2,
                ),
            )
                if field_name in output_variables
                    field = getfield(integrals, field_name)[ii, jj, kk]

                    HDF5.set_extent_dims(
                        file[string(field_name)],
                        (x_size, y_size, z_size, iout),
                    )
                    file[string(field_name)][iid, jjd, kkd, iout] =
                        @share field[ii, jj, kk] * scaling
                end
            end

            # Write GW tendencies.
            for (field_name, scaling) in zip(
                (:dudt, :dvdt, :dthetadt),
                (
                    rhoref * uref / tref,
                    rhoref * uref / tref,
                    rhoref * thetaref / tref,
                ),
            )
                if field_name in output_variables
                    field = getfield(tendencies, field_name)[ii, jj, kk]

                    HDF5.set_extent_dims(
                        file[string(field_name)],
                        (x_size, y_size, z_size, iout),
                    )
                    file[string(field_name)][iid, jjd, kkd, iout] =
                        @share field[ii, jj, kk] * scaling
                end
            end

            # Write elastic-mode-selection data.
            if elastic_mode_selection && ko == 0
                for field in (:launch_mode_count, :launch_power_fraction)
                    if field in output_variables
                        HDF5.set_extent_dims(
                            file[string(field)],
                            (x_size, y_size, iout),
                        )
                        file[string(field)][iid, jjd, iout] =
                            getfield(state.wkb.elastic_mode_selection, field)[
                                ii,
                                jj,
                            ]
                    end
                end
            end
        end

        return
    end

    return iout
end
