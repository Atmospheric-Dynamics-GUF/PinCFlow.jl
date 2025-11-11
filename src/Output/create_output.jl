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
        attributes(dset)["long_name"] = "x-coordinates"

        dset =
            create_dataset(file, "y", datatype(Float32), dataspace((y_size,)))
        attributes(dset)["units"] = "m"
        attributes(dset)["long_name"] = "y-coordinates"

        dset = create_dataset(
            file,
            "z",
            datatype(Float32),
            dataspace((x_size, y_size, z_size));
            chunk = (cx, cy, cz),
        )
        attributes(dset)["units"] = "m"
        attributes(dset)["long_name"] = "z-coordinates"

        dset = create_dataset(
            file,
            "ztilde",
            datatype(Float32),
            dataspace((x_size, y_size, z_size + 1));
            chunk = (cx, cy, cz),
        )
        attributes(dset)["units"] = "m"
        attributes(dset)["long_name"] = "z-coordinates"

        dset = create_dataset(
            file,
            "t",
            datatype(Float32),
            dataspace((0,), (-1,));
            chunk = (ct,),
        )
        attributes(dset)["units"] = "s"
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
            attributes(dset)["units"] = "kg m^-3"
            attributes(dset)["long_name"] = "background density"

            dset = create_dataset(
                file,
                "thetabar",
                datatype(Float32),
                dataspace((x_size, y_size, z_size));
                chunk = (cx, cy, cz),
            )
            attributes(dset)["units"] = "K"
            attributes(dset)["long_name"] = "background potential temperature"

            dset = create_dataset(
                file,
                "n2",
                datatype(Float32),
                dataspace((x_size, y_size, z_size));
                chunk = (cx, cy, cz),
            )
            attributes(dset)["units"] = "s^-2"
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
                attributes(dset)["units"] = "kg K m^-3"
                attributes(dset)["long_name"] = "mass-weighted potential temperature"
            else
                dset = create_dataset(
                    file,
                    "p",
                    datatype(Float32),
                    dataspace((x_size, y_size, z_size));
                    chunk = (cx, cy, cz),
                )
                attributes(dset)["units"] = "kg K m^-3"
                attributes(dset)["long_name"] = "mass-weighted potential temperature"
            end
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
            attributes(dset)["units"] = "kg m^-3"
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
            attributes(dset)["units"] = "m s^-1"
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
            attributes(dset)["units"] = "m s^-1"
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
            attributes(dset)["units"] = "m s^-1"
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
            attributes(dset)["units"] = "m s^-1"
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
            attributes(dset)["units"] = "m s^-1"
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
            attributes(dset)["units"] = "m s^-1"
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
            attributes(dset)["units"] = "m s^-1"
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
            attributes(dset)["units"] = "m s^-1"
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
            attributes(dset)["long_name"] = "potential temperature fluctuations"
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
            attributes(dset)["units"] = ""
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
                attributes(dset)["units"] = ""
                attributes(dset)["long_name"] = "tracer mixing ratio"
            end

            if state.namelists.tracer.leading_order_impact &&
               :dchidt in output_variables
                for (ifld, field) in enumerate(fieldnames(TracerWKBImpact))
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
                    attributes(dset)["units"] =
                        ["m s^-1", "m s^-1", "m s^-1", "s^-1"][ifld]
                    attributes(dset)["long_name"] = [
                        "zonal GW-tracer flux",
                        "meridiona GW-tracer flux",
                        "vertical GW-tracer flux",
                        "gravity-wave-tracer flux convergence",
                    ][ifld]
                end
            end
        end

        # Create datasets for WKB variables.
        if wkb_mode != NoWKB()

            # Create datasets for ray-volume properties.
            if prepare_restart || save_ray_volumes
                for (ifld, field) in enumerate((
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
                ))
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
                    attributes(dset)["units"] = [
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
                        "kg s^-1",
                    ][ifld]
                    attributes(dset)["long_name"] = [
                        "ray volume position in x-direction",
                        "ray volume position in y-direction",
                        "ray volume position in z-direction",
                        "ray volume extent in x-direction",
                        "ray volume extent in y-direction",
                        "ray volume extent in z-direction",
                        "ray volume position in k-direction",
                        "ray volume position in l-direction",
                        "ray volume position in m-direction",
                        "ray volume extent in k-direction",
                        "ray volume extent in l-direction",
                        "ray volume extent in m-direction",
                        "ray volume wave-action density",
                    ][ifld]
                end
            end

            # Create datasets for GW tendencies.
            for (ifld, field) in enumerate((:dudt, :dvdt, :dthetadt))
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
                    attributes(dset)["units"] =
                        ["m s^-2", "m s^-2", "K s^-1"][ifld]
                    attributes(dset)["long_name"] = [
                        "zonal wind GW tendency",
                        "meridional wind GW tendency",
                        "potential temperature GW tendency",
                    ][ifld]
                end
            end
        end

        return
    end

    return
end
