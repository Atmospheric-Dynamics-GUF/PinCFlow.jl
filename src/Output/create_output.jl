"""
```julia
create_output(state::State)
```

Create an HDF5 output file with one dataset for each variable.

The dimensions of the datasets are set to those of the domain, whereas the chunks are set to the dimensions of the MPI subdomains, in preparation for parallel output. Datasets for the grid, i.e. the fields `x`, `y` and `zc` of `state.grid`, the time and the fields of `state.atmosphere` are always created, regardless of the specifications in `state.namelists.output`. The one exception to this is the Boussinesq mode, in which no datasets are created for the fields of `state.atmosphere`, since they do not have a spatial dependence.

# Arguments

  - `state`: Model state.
"""
function create_output end

function create_output(state::State, machine_start_time::DateTime)
    (; x_size, y_size, z_size, npx, npy, npz) = state.namelists.domain
    (; prepare_restart, save_ray_volumes, output_variables, output_file) =
        state.namelists.output
    (; model) = state.namelists.atmosphere
    (; wkb_mode) = state.namelists.wkb
    (; comm) = state.domain
    (; nray_max) = state.wkb

    # Set the chunk dimensions.
    cr = nray_max
    cx = div(x_size, npx)
    cy = div(y_size, npy)
    cz = div(z_size, npz)
    ct = 1

    # Prepare the output.
    h5open(output_file, "w", comm) do file
        attributes(file)["Title"] = "PinCFlow.jl data"
        attributes(
            file,
        )["Institution"] = "Institute for Atmospheric and Environmental Sciences, Goethe University Frankfurt, Germany"
        attributes(file)["Date"] = string(Dates.Date(machine_start_time))
        attributes(file)["Time"] = string(Dates.Time(machine_start_time))

        # Create datasets for the dimensions.
        dset =
            create_dataset(file, "x", datatype(Float32), dataspace((x_size,)))
        attributes(dset)["units"] = "m"
        attributes(dset)["label"] = L"x\,[\mathrm{m}]"
        attributes(dset)["long_name"] = "x-coordinates"

        dset =
            create_dataset(file, "y", datatype(Float32), dataspace((y_size,)))
        attributes(dset)["units"] = "m"
        attributes(dset)["label"] = L"y\,[\mathrm{m}]"
        attributes(dset)["long_name"] = "y-coordinates"

        dset = create_dataset(
            file,
            "z",
            datatype(Float32),
            dataspace((x_size, y_size, z_size));
            chunk = (cx, cy, cz),
        )
        attributes(dset)["units"] = "m"
        attributes(dset)["label"] = L"z\,[\mathrm{m}]"
        attributes(dset)["long_name"] = "z-coordinates"

        dset = create_dataset(
            file,
            "ztilde",
            datatype(Float32),
            dataspace((x_size, y_size, z_size + 1));
            chunk = (cx, cy, cz),
        )
        attributes(dset)["units"] = "m"
        attributes(dset)["label"] = L"z_{\mathrm{s}}\,[\mathrm{m}]"
        attributes(dset)["long_name"] = "staggered z-coordinates"

        dset = create_dataset(
            file,
            "t",
            datatype(Float32),
            dataspace((0,), (-1,));
            chunk = (ct,),
        )
        attributes(dset)["units"] = "s"
        attributes(dset)["label"] = L"t\,[\mathrm{s}]"
        attributes(dset)["long_name"] = "time"

        # Create datasets for the background.
        if model != Boussinesq()
            dset = create_dataset(
                file,
                "rhobar",
                datatype(Float32),
                dataspace((x_size, y_size, z_size));
                chunk = (cx, cy, cz),
            )
            attributes(dset)["units"] = "kg*m^-3"
            attributes(
                dset,
            )["label"] = L"\overline{\rho}\,[\mathrm{kg\,m^{-3}}]"
            attributes(dset)["long_name"] = "density background"

            dset = create_dataset(
                file,
                "thetabar",
                datatype(Float32),
                dataspace((x_size, y_size, z_size));
                chunk = (cx, cy, cz),
            )
            attributes(dset)["units"] = "K"
            attributes(dset)["label"] = L"\overline{\theta}\,[\mathrm{K}]"
            attributes(dset)["long_name"] = "potential-temperature background"

            dset = create_dataset(
                file,
                "n2",
                datatype(Float32),
                dataspace((x_size, y_size, z_size));
                chunk = (cx, cy, cz),
            )
            attributes(dset)["units"] = "s^-2"
            attributes(dset)["label"] = L"N^2\,[\mathrm{s^{-2}}]"
            attributes(dset)["long_name"] = "squared buoyancy frequency"

            if model == Compressible()
                dset = create_dataset(
                    file,
                    "p",
                    datatype(Float32),
                    dataspace(
                        (x_size, y_size, z_size, 0),
                        (x_size, y_size, z_size, -1),
                    );
                    chunk = (cx, cy, cz, ct),
                )
            else
                dset = create_dataset(
                    file,
                    "p",
                    datatype(Float32),
                    dataspace((x_size, y_size, z_size));
                    chunk = (cx, cy, cz),
                )
            end
            attributes(dset)["units"] = "kg*K*m^-3"
            attributes(dset)["label"] = L"P\,[\mathrm{kg\,K\,m^{-3}}]"
            attributes(dset)["long_name"] = "mass-weighted potential temperature"
        end

        # Create datasets for the prognostic variables.
        if prepare_restart || :rhop in output_variables
            dset = create_dataset(
                file,
                "rhop",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "kg*m^-3"
            attributes(dset)["label"] = L"\rho'\,[\mathrm{kg\,m^{-3}}]"
            attributes(dset)["long_name"] = "density fluctuations"
        end
        if :u in output_variables
            dset = create_dataset(
                file,
                "u",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "m*s^-1"
            attributes(dset)["label"] = L"u\,[\mathrm{m\,s^{-1}}]"
            attributes(dset)["long_name"] = "zonal wind"
        end
        if prepare_restart || :us in output_variables
            dset = create_dataset(
                file,
                "us",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "m*s^-1"
            attributes(dset)["label"] = L"u_{\mathrm{s}}\,[\mathrm{m\,s^{-1}}]"
            attributes(dset)["long_name"] = "staggered zonal wind"
        end
        if :v in output_variables
            dset = create_dataset(
                file,
                "v",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "m*s^-1"
            attributes(dset)["label"] = L"v\,[\mathrm{m\,s^{-1}}]"
            attributes(dset)["long_name"] = "meridional wind"
        end
        if prepare_restart || :vs in output_variables
            dset = create_dataset(
                file,
                "vs",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "m*s^-1"
            attributes(dset)["label"] = L"v_{\mathrm{s}}\,[\mathrm{m\,s^{-1}}]"
            attributes(dset)["long_name"] = "staggered meridional wind"
        end
        if :w in output_variables
            dset = create_dataset(
                file,
                "w",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "m*s^-1"
            attributes(dset)["label"] = L"w\,[\mathrm{m\,s^{-1}}]"
            attributes(dset)["long_name"] = "vertical wind"
        end
        if :ws in output_variables
            dset = create_dataset(
                file,
                "ws",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "m*s^-1"
            attributes(dset)["label"] = L"w_{\mathrm{s}}\,[\mathrm{m\,s^{-1}}]"
            attributes(dset)["long_name"] = "staggered vertical wind"
        end
        if :wt in output_variables
            dset = create_dataset(
                file,
                "wt",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "m*s^-1"
            attributes(dset)["label"] = L"\widehat{w}\,[\mathrm{m\,s^{-1}}]"
            attributes(dset)["long_name"] = "transformed vertical wind"
        end
        if prepare_restart || :wts in output_variables
            dset = create_dataset(
                file,
                "wts",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "m*s^-1"
            attributes(
                dset,
            )["label"] = L"\widehat{w}_{\mathrm{s}}\,[\mathrm{m\,s^{-1}}]"
            attributes(dset)["long_name"] = "staggered transformed vertical wind"
        end
        if :thetap in output_variables
            dset = create_dataset(
                file,
                "thetap",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "K"
            attributes(dset)["label"] = L"\theta'\,[\mathrm{K}]"
            attributes(dset)["long_name"] = "potential-temperature fluctuations"
        end
        if prepare_restart || :pip in output_variables
            dset = create_dataset(
                file,
                "pip",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
            attributes(dset)["units"] = "1"
            attributes(dset)["label"] = L"\pi'"
            attributes(dset)["long_name"] = "Exner-pressure fluctuations"
        end

        if !(typeof(state.namelists.tracer.tracer_setup) <: NoTracer)
            for field in fieldnames(TracerPredictands)
                dset = create_dataset(
                    file,
                    string(field),
                    datatype(Float32),
                    dataspace(
                        (x_size, y_size, z_size, 0),
                        (x_size, y_size, z_size, -1),
                    );
                    chunk = (cx, cy, cz, ct),
                )
                attributes(dset)["units"] = "1"
                attributes(dset)["label"] = L"\chi"
                attributes(dset)["long_name"] = "tracer mixing ratio"
            end

            if state.namelists.tracer.leading_order_impact &&
               :dchidt in output_variables
                for (field, units, label, long_name) in zip(
                    fieldnames(TracerWKBImpact),
                    ("m*s^-1", "m*s^-1", "m*s^-1", "s^-1"),
                    (
                        L"\langle u'\chi' \rangle\,[\mathrm{m\,s^{-1}}]",
                        L"\langle v'\chi' \rangle\,[\mathrm{m\,s^{-1}}]",
                        L"\langle w'\chi' \rangle\,[\mathrm{m\,s^{-1}}]",
                        L"(\partial_t \chi_\mathrm{b})^{(0)}_\mathrm{w},[\mathrm{s^{-1}}]",
                    ),
                    (
                        "zonal GW-tracer flux",
                        "meridional GW-tracer flux",
                        "vertical GW-tracer flux",
                        "GW-tracer flux convergence",
                    ),
                )
                    dset = create_dataset(
                        file,
                        string(field),
                        datatype(Float32),
                        dataspace(
                            (x_size, y_size, z_size, 0),
                            (x_size, y_size, z_size, -1),
                        );
                        chunk = (cx, cy, cz, ct),
                    )
                    attributes(dset)["units"] = units
                    attributes(dset)["label"] = label
                    attributes(dset)["long_name"] = long_name
                end
            end
        end

        # Create datasets for WKB variables.
        if wkb_mode != NoWKB()

            # Create datasets for ray-volume properties.
            if prepare_restart || save_ray_volumes
                if x_size == 1 && y_size == 1
                    nr_units = "kg*s^-1"
                    nr_label = L"\mathcal{N}_r\,[\mathrm{kg\,s^{-1}}]"
                elseif x_size > 1 && y_size > 1
                    nr_units = "kg*m^2*s^-1"
                    nr_label = L"\mathcal{N}_r\,[\mathrm{kg\,m^2\,s^{-1}}]"
                else
                    nr_units = "kg*m*s^-1"
                    nr_label = L"\mathcal{N}_r\,[\mathrm{kg\,m\,s^{-1}}]"
                end
                for (field, units, label, long_name) in zip(
                    (
                        "xr",
                        "yr",
                        "zr",
                        "dxr",
                        "dyr",
                        "dzr",
                        "kr",
                        "lr",
                        "mr",
                        "dkr",
                        "dlr",
                        "dmr",
                        "nr",
                    ),
                    (
                        "m",
                        "m",
                        "m",
                        "m",
                        "m",
                        "m",
                        "m^-1",
                        "m^-1",
                        "m^-1",
                        "m^-1",
                        "m^-1",
                        "m^-1",
                        nr_units,
                    ),
                    (
                        L"x_{r}\,[\mathrm{m}]",
                        L"y_{r}\,[\mathrm{m}]",
                        L"z_{r}\,[\mathrm{m}]",
                        L"\Delta x_{r}\,[\mathrm{m}]",
                        L"\Delta y_{r}\,[\mathrm{m}]",
                        L"\Delta z_{r}\,[\mathrm{m}]",
                        L"k_{r}\,[\mathrm{m^{-1}}]",
                        L"l_{r}\,[\mathrm{m^{-1}}]",
                        L"m_{r}\,[\mathrm{m^{-1}}]",
                        L"\Delta k_{r}\,[\mathrm{m^{-1}}]",
                        L"\Delta l_{r}\,[\mathrm{m^{-1}}]",
                        L"\Delta m_{r}\,[\mathrm{m^{-1}}]",
                        nr_label,
                    ),
                    (
                        "ray-volume position in x",
                        "ray-volume position in y",
                        "ray-volume position in z",
                        "ray-volume extent in x",
                        "ray-volume extent in y",
                        "ray-volume extent in z",
                        "ray-volume position in k",
                        "ray-volume position in l",
                        "ray-volume position in m",
                        "ray-volume extent in k",
                        "ray-volume extent in l",
                        "ray-volume extent in m",
                        "ray-volume phase-space wave-action density",
                    ),
                )
                    dset = create_dataset(
                        file,
                        field,
                        datatype(Float32),
                        dataspace(
                            (nray_max, x_size, y_size, z_size + 1, 0),
                            (nray_max, x_size, y_size, z_size + 1, -1),
                        );
                        chunk = (cr, cx, cy, cz, ct),
                    )
                    attributes(dset)["units"] = units
                    attributes(dset)["label"] = label
                    attributes(dset)["long_name"] = long_name
                end
            end

            # Create datasets for GW tendencies.
            for (field, units, label, long_name) in zip(
                (:dudt, :dvdt, :dthetadt),
                ("kg*m^-2*s^-2", "kg*m^-2*s^-2", "kg*K*m^-3*s^-1"),
                (
                    L"[\partial_t (\rho_\mathrm{b} u_\mathrm{b})]_\mathrm{w}\,[\mathrm{kg\,m^{-2}\,s^{-2}}]",
                    L"[\partial_t (\rho_\mathrm{b} v_\mathrm{b})]_\mathrm{w}\,[\mathrm{kg\,m^{-2}\,s^{-2}}]",
                    L"[\partial_t (P_\mathrm{b})]_\mathrm{w}\,[\mathrm{kg\,K\,m^{-3}\,s^{-1}}]",
                ),
                (
                    "zonal-momentum GW forcing",
                    "meridional-momentum GW forcing",
                    "mass-weighted potential-temperature GW forcing",
                ),
            )
                if field in output_variables
                    dset = create_dataset(
                        file,
                        string(field),
                        datatype(Float32),
                        dataspace(
                            (x_size, y_size, z_size, 0),
                            (x_size, y_size, z_size, -1),
                        );
                        chunk = (cx, cy, cz, ct),
                    )
                    attributes(dset)["units"] = units
                    attributes(dset)["label"] = label
                    attributes(dset)["long_name"] = long_name
                end
            end
        end

        return
    end

    return
end
