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
    (; comm, master) = state.domain
    (; nray_max) = state.wkb

    # Set the chunk dimensions.
    cr = nray_max
    cx = div(x_size, npx)
    cy = div(y_size, npy)
    cz = div(z_size, npz)
    ct = 1

    # Create the output file and the datasets.
    h5open(output_file, "w", comm) do file

        # Create datasets for the dimensions.
        create_dataset(file, "x", datatype(Float32), dataspace((x_size,)))
        create_dataset(file, "y", datatype(Float32), dataspace((y_size,)))
        create_dataset(
            file,
            "z",
            datatype(Float32),
            dataspace((x_size, y_size, z_size));
            chunk = (cx, cy, cz),
        )
        create_dataset(
            file,
            "ztilde",
            datatype(Float32),
            dataspace((x_size, y_size, z_size + 1));
            chunk = (cx, cy, cz),
        )
        create_dataset(
            file,
            "t",
            datatype(Float32),
            dataspace((0,), (-1,));
            chunk = (ct,),
        )

        # Create datasets for the background.
        if model != Boussinesq()
            create_dataset(
                file,
                "rhobar",
                datatype(Float32),
                dataspace((x_size, y_size, z_size));
                chunk = (cx, cy, cz),
            )
            create_dataset(
                file,
                "thetabar",
                datatype(Float32),
                dataspace((x_size, y_size, z_size));
                chunk = (cx, cy, cz),
            )
            create_dataset(
                file,
                "n2",
                datatype(Float32),
                dataspace((x_size, y_size, z_size));
                chunk = (cx, cy, cz),
            )
            if model == Compressible()
                create_dataset(
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
                create_dataset(
                    file,
                    "p",
                    datatype(Float32),
                    dataspace((x_size, y_size, z_size));
                    chunk = (cx, cy, cz),
                )
            end
        end

        # Create datasets for the prognostic variables.

        if prepare_restart || :rhop in output_variables
            create_dataset(
                file,
                "rhop",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if :u in output_variables
            create_dataset(
                file,
                "u",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if prepare_restart || :us in output_variables
            create_dataset(
                file,
                "us",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if :v in output_variables
            create_dataset(
                file,
                "v",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if prepare_restart || :vs in output_variables
            create_dataset(
                file,
                "vs",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if :w in output_variables
            create_dataset(
                file,
                "w",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if :ws in output_variables
            create_dataset(
                file,
                "ws",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if :wt in output_variables
            create_dataset(
                file,
                "wt",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if prepare_restart || :wts in output_variables
            create_dataset(
                file,
                "wts",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if :thetap in output_variables
            create_dataset(
                file,
                "thetap",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if prepare_restart || :pip in output_variables
            create_dataset(
                file,
                "pip",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if !(typeof(state.namelists.tracer.tracer_setup) <: NoTracer)
            for field in fieldnames(TracerPredictands)
                create_dataset(
                    file,
                    string(field),
                    datatype(Float32),
                    dataspace(
                        (x_size, y_size, z_size, 0),
                        (x_size, y_size, z_size, -1),
                    );
                    chunk = (cx, cy, cz, ct),
                )
            end

            if state.namelists.tracer.leading_order_impact &&
               :dchidt in output_variables
                for field in fieldnames(TracerWKBImpact)
                    create_dataset(
                        file,
                        string(field),
                        datatype(Float32),
                        dataspace(
                            (x_size, y_size, z_size, 0),
                            (x_size, y_size, z_size, -1),
                        );
                        chunk = (cx, cy, cz, ct),
                    )
                end
            end
        end

        # Create datasets for WKB variables.
        if wkb_mode != NoWKB()

            # Create datasets for ray-volume properties.
            if prepare_restart || save_ray_volumes
                for field in (
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
                )
                    create_dataset(
                        file,
                        field,
                        datatype(Float32),
                        dataspace(
                            (nray_max, x_size, y_size, z_size + 1, 0),
                            (nray_max, x_size, y_size, z_size + 1, -1),
                        );
                        chunk = (cr, cx, cy, cz, ct),
                    )
                end
            end

            # Create datasets for GW tendencies.
            for field in (:dudt, :dvdt, :dthetadt)
                if field in output_variables
                    create_dataset(
                        file,
                        string(field),
                        datatype(Float32),
                        dataspace(
                            (x_size, y_size, z_size, 0),
                            (x_size, y_size, z_size, -1),
                        );
                        chunk = (cx, cy, cz, ct),
                    )
                end
            end
        end

        return
    end

    # Add attributes and namelists.
    master && h5open(output_file, "r+") do file
        attributes(file)["Title"] = "PinCFlow.jl data"
        attributes(
            file,
        )["Institution"] = "Institute for Atmospheric and Environmental Sciences, Goethe University Frankfurt, Germany"
        attributes(file)["Date"] = string(Dates.Date(machine_start_time))
        attributes(file)["Time"] = string(Dates.Time(machine_start_time))

        attributes(file["x"])["units"] = "m"
        attributes(file["x"])["label"] = L"x\ [\mathrm{m}]"
        attributes(file["x"])["long_name"] = "x-coordinates"

        attributes(file["y"])["units"] = "m"
        attributes(file["y"])["label"] = L"y\ [\mathrm{m}]"
        attributes(file["y"])["long_name"] = "y-coordinates"

        attributes(file["z"])["units"] = "m"
        attributes(file["z"])["label"] = L"z\ [\mathrm{m}]"
        attributes(file["z"])["long_name"] = "z-coordinates"

        attributes(file["ztilde"])["units"] = "m"
        attributes(file["ztilde"])["label"] =
            L"z_{\mathrm{s}}\ [\mathrm{m}]"
        attributes(file["ztilde"])["long_name"] = "staggered z-coordinates"

        attributes(file["t"])["units"] = "s"
        attributes(file["t"])["label"] = L"t\ [\mathrm{s}]"
        attributes(file["t"])["long_name"] = "time"

        if model != Boussinesq()
            attributes(file["rhobar"])["units"] = "kg*m^-3"
            attributes(file["rhobar"])["label"] =
                L"\overline{\rho}\ [\mathrm{kg\ m^{-3}}]"
            attributes(file["rhobar"])["long_name"] = "density background"

            attributes(file["thetabar"])["units"] = "K"
            attributes(file["thetabar"])["label"] =
                L"\overline{\theta}\ [\mathrm{K}]"
            attributes(file["thetabar"])["long_name"] = "potential-temperature background"

            attributes(file["n2"])["units"] = "s^-2"
            attributes(file["n2"])["label"] = L"N^2\ [\mathrm{s^{-2}}]"
            attributes(file["n2"])["long_name"] = "squared buoyancy frequency"

            attributes(file["p"])["units"] = "kg*K*m^-3"
            attributes(file["p"])["label"] = L"P\ [\mathrm{kg\ K\ m^{-3}}]"
            attributes(file["p"])["long_name"] = "mass-weighted potential temperature"
        end

        if prepare_restart || :rhop in output_variables
            attributes(file["rhop"])["units"] = "kg*m^-3"
            attributes(file["rhop"])["label"] =
                L"\rho'\ [\mathrm{kg\ m^{-3}}]"
            attributes(file["rhop"])["long_name"] = "density fluctuations"
        end

        if :u in output_variables
            attributes(file["u"])["units"] = "m*s^-1"
            attributes(file["u"])["label"] = L"u\ [\mathrm{m\ s^{-1}}]"
            attributes(file["u"])["long_name"] = "zonal wind"
        end

        if prepare_restart || :us in output_variables
            attributes(file["us"])["units"] = "m*s^-1"
            attributes(file["us"])["label"] =
                L"u_{\mathrm{s}}\ [\mathrm{m\ s^{-1}}]"
            attributes(file["us"])["long_name"] = "staggered zonal wind"
        end

        if :v in output_variables
            attributes(file["v"])["units"] = "m*s^-1"
            attributes(file["v"])["label"] = L"v\ [\mathrm{m\ s^{-1}}]"
            attributes(file["v"])["long_name"] = "meridional wind"
        end

        if prepare_restart || :vs in output_variables
            attributes(file["vs"])["units"] = "m*s^-1"
            attributes(file["vs"])["label"] =
                L"v_{\mathrm{s}}\ [\mathrm{m\ s^{-1}}]"
            attributes(file["vs"])["long_name"] = "staggered meridional wind"
        end

        if :w in output_variables
            attributes(file["w"])["units"] = "m*s^-1"
            attributes(file["w"])["label"] = L"w\ [\mathrm{m\ s^{-1}}]"
            attributes(file["w"])["long_name"] = "vertical wind"
        end

        if :ws in output_variables
            attributes(file["ws"])["units"] = "m*s^-1"
            attributes(file["ws"])["label"] =
                L"w_{\mathrm{s}}\ [\mathrm{m\ s^{-1}}]"
            attributes(file["ws"])["long_name"] = "staggered vertical wind"
        end

        if :wt in output_variables
            attributes(file["wt"])["units"] = "m*s^-1"
            attributes(file["wt"])["label"] =
                L"\widehat{w}\ [\mathrm{m\ s^{-1}}]"
            attributes(file["wt"])["long_name"] = "transformed vertical wind"
        end

        if prepare_restart || :wts in output_variables
            attributes(file["wts"])["units"] = "m*s^-1"
            attributes(file["wts"])["label"] =
                L"\widehat{w}_{\mathrm{s}}\ [\mathrm{m\ s^{-1}}]"
            attributes(file["wts"])["long_name"] = "staggered transformed vertical wind"
        end

        if :thetap in output_variables
            attributes(file["thetap"])["units"] = "K"
            attributes(file["thetap"])["label"] = L"\theta'\ [\mathrm{K}]"
            attributes(file["thetap"])["long_name"] = "potential-temperature fluctuations"
        end

        if prepare_restart || :pip in output_variables
            attributes(file["pip"])["units"] = "1"
            attributes(file["pip"])["label"] = L"\pi'"
            attributes(file["pip"])["long_name"] = "Exner-pressure fluctuations"
        end

        if !(typeof(state.namelists.tracer.tracer_setup) <: NoTracer)
            for field in fieldnames(TracerPredictands)
                attributes(file[string(field)])["units"] = "1"
                attributes(file[string(field)])["label"] = L"\chi"
                attributes(file[string(field)])["long_name"] = "tracer mixing ratio"
            end

            if state.namelists.tracer.leading_order_impact &&
               :dchidt in output_variables
                for (field, units, label, long_name) in zip(
                    fieldnames(TracerWKBImpact),
                    ("m*s^-1", "m*s^-1", "m*s^-1", "s^-1"),
                    (
                        L"\langle u'\chi' \rangle\ [\mathrm{m\ s^{-1}}]",
                        L"\langle v'\chi' \rangle\ [\mathrm{m\ s^{-1}}]",
                        L"\langle w'\chi' \rangle\ [\mathrm{m\ s^{-1}}]",
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

        if !(
            typeof(state.namelists.turbulence.turbulence_scheme) <:
            NoTurbulence
        )
            for field in fieldnames(TurbulencePredictands)
                create_dataset(
                    file,
                    string(field),
                    datatype(Float32),
                    dataspace(
                        (x_size, y_size, z_size, 0),
                        (x_size, y_size, z_size, -1),
                    );
                    chunk = (cx, cy, cz, ct),
                )
            end

            create_dataset(
                file,
                "shear-production",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )

            create_dataset(
                file,
                "buoyancy-production",
                datatype(Float32),
                dataspace(
                    (x_size, y_size, z_size, 0),
                    (x_size, y_size, z_size, -1),
                );
                chunk = (cx, cy, cz, ct),
            )
        end

        if wkb_mode != NoWKB()
            if prepare_restart || save_ray_volumes
                if x_size == 1 && y_size == 1
                    nr_units = "kg*s^-1"
                    nr_label = L"\mathcal{N}_r\ [\mathrm{kg\ s^{-1}}]"
                elseif x_size > 1 && y_size > 1
                    nr_units = "kg*m^2*s^-1"
                    nr_label = L"\mathcal{N}_r\ [\mathrm{kg\ m^2\ s^{-1}}]"
                else
                    nr_units = "kg*m*s^-1"
                    nr_label = L"\mathcal{N}_r\ [\mathrm{kg\ m\ s^{-1}}]"
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
                        L"x_{r}\ [\mathrm{m}]",
                        L"y_{r}\ [\mathrm{m}]",
                        L"z_{r}\ [\mathrm{m}]",
                        L"\Delta x_{r}\ [\mathrm{m}]",
                        L"\Delta y_{r}\ [\mathrm{m}]",
                        L"\Delta z_{r}\ [\mathrm{m}]",
                        L"k_{r}\ [\mathrm{m^{-1}}]",
                        L"l_{r}\ [\mathrm{m^{-1}}]",
                        L"m_{r}\ [\mathrm{m^{-1}}]",
                        L"\Delta k_{r}\ [\mathrm{m^{-1}}]",
                        L"\Delta l_{r}\ [\mathrm{m^{-1}}]",
                        L"\Delta m_{r}\ [\mathrm{m^{-1}}]",
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
                    attributes(file[field])["units"] = units
                    attributes(file[field])["label"] = label
                    attributes(file[field])["long_name"] = long_name
                end
            end

            # Create datasets for GW tendencies.
            for (field, units, label, long_name) in zip(
                (:dudt, :dvdt, :dthetadt),
                ("kg*m^-2*s^-2", "kg*m^-2*s^-2", "kg*K*m^-3*s^-1"),
                (
                    L"[\partial_t (\rho_\mathrm{b} u_\mathrm{b})]_\mathrm{w}\ [\mathrm{kg\ m^{-2}\ s^{-2}}]",
                    L"[\partial_t (\rho_\mathrm{b} v_\mathrm{b})]_\mathrm{w}\ [\mathrm{kg\ m^{-2}\ s^{-2}}]",
                    L"[\partial_t (P_\mathrm{b})]_\mathrm{w}\ [\mathrm{kg\ K\ m^{-3}\ s^{-1}}]",
                ),
                (
                    "zonal-momentum GW forcing",
                    "meridional-momentum GW forcing",
                    "mass-weighted potential-temperature GW forcing",
                ),
            )
                if field in output_variables
                    attributes(file[string(field)])["units"] = units
                    attributes(file[string(field)])["label"] = label
                    attributes(file[string(field)])["long_name"] = long_name
                end
            end
        end

        create_group(file, "namelists")
        for namelist in fieldnames(Namelists)
            create_group(file["namelists"], string(namelist))
            for parameter in
                fieldnames(typeof(getfield(state.namelists, namelist)))
                value = getfield(getfield(state.namelists, namelist), parameter)
                file["namelists"][string(namelist)][string(parameter)] =
                    typeof(value) <: AbstractString ? "\"" * value * "\"" :
                    string(value)
            end
        end

        return
    end

    MPI.Barrier(comm)

    return
end
