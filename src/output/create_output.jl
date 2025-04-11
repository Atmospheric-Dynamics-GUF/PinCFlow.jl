function create_output(state::State)
    # Get all necessary fields.
    (; sizex, sizey, sizez) = state.namelists.domain
    (; prepare_restart, output_variables, folder) = state.namelists.output
    (; comm, nx, ny, nz) = state.domain

    # Prepare the output.
    h5open(folder * "/pincflow_output.h5", "w", comm) do file

        # Create datasets for the dimensions.
        create_dataset(file, "x", datatype(Float32), dataspace((sizex,)))
        create_dataset(file, "y", datatype(Float32), dataspace((sizey,)))
        create_dataset(
            file,
            "z",
            datatype(Float32),
            dataspace((sizex, sizey, sizez));
            chunk = (nx, ny, nz),
        )
        create_dataset(
            file,
            "t",
            datatype(Float32),
            dataspace((0,), (-1,));
            chunk = (1,),
        )

        # Create datasets for the background.
        for label in ("rhobar", "thetabar", "n2", "p")
            create_dataset(
                file,
                label,
                datatype(Float32),
                dataspace((sizex, sizey, sizez));
                chunk = (nx, ny, nz),
            )
        end

        # Create datasets for the prognostic variables.
        if prepare_restart || :rhop in output_variables
            create_dataset(
                file,
                "rhop",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if :u in output_variables
            create_dataset(
                file,
                "u",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if prepare_restart || :us in output_variables
            create_dataset(
                file,
                "us",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if :v in output_variables
            create_dataset(
                file,
                "v",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if prepare_restart || :vs in output_variables
            create_dataset(
                file,
                "vs",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if :w in output_variables
            create_dataset(
                file,
                "w",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if :ws in output_variables
            create_dataset(
                file,
                "ws",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if :wtfc in output_variables
            create_dataset(
                file,
                "wtfc",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if prepare_restart || :wstfc in output_variables
            create_dataset(
                file,
                "wstfc",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if :thetap in output_variables
            create_dataset(
                file,
                "thetap",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end
        if prepare_restart || :pip in output_variables
            create_dataset(
                file,
                "pip",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (nx, ny, nz, 1),
            )
        end

        # Return.
        return
    end

    # Return.
    return
end
