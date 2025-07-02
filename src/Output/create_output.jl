"""
    create_output(state::State)

Initialize HDF5 output file with dimensions, attributes, and datasets.

Creates output file structure with grid coordinates, variable datasets, and
metadata. Sets up chunking and compression for efficient parallel I/O.

# Arguments

  - `state::State`: Simulation state containing grid and configuration

# File Structure

  - **Dimensions**: `nx`, `ny`, `nz`, `nt` for spatial and temporal extents
  - **Coordinates**: `x`, `y`, `z`, `t` arrays with proper scaling
  - **Variables**: Datasets for all prognostic and diagnostic fields
  - **Attributes**: Simulation metadata, constants, and namelist parameters

# Implementation

  - **Chunking**: Optimized for domain decomposition access patterns
  - **Compression**: Applies deflate compression for storage efficiency
  - **Parallel**: Creates collective datasets for MPI-parallel writing
"""
function create_output(state::State)
    (; sizex, sizey, sizez, npx, npy, npz) = state.namelists.domain
    (; prepare_restart, save_ray_volumes, output_variables, output_file) =
        state.namelists.output
    (; model, testcase) = state.namelists.setting
    (; comm) = state.domain
    (; nray_max) = state.wkb

    # Set the chunk dimensions.
    cr = nray_max
    cx = div(sizex, npx)
    cy = div(sizey, npy)
    cz = div(sizez, npz)
    ct = 1

    # Prepare the output.
    h5open(output_file, "w", comm) do file

        # Create datasets for the dimensions.
        create_dataset(file, "x", datatype(Float32), dataspace((sizex,)))
        create_dataset(file, "y", datatype(Float32), dataspace((sizey,)))
        create_dataset(
            file,
            "z",
            datatype(Float32),
            dataspace((sizex, sizey, sizez));
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
            for label in ("rhobar", "thetabar")
                create_dataset(
                    file,
                    label,
                    datatype(Float32),
                    dataspace((sizex, sizey, sizez));
                    chunk = (cx, cy, cz),
                )
            end

            if model == Compressible()
                for label in ("n2", "p")
                    create_dataset(
                        file,
                        label,
                        datatype(Float32),
                        dataspace(
                            (sizex, sizey, sizez, 0),
                            (sizex, sizey, sizez, -1),
                        );
                        chunk = (cx, cy, cz, ct),
                    )
                end
            else
                for label in ("n2", "p")
                    create_dataset(
                        file,
                        label,
                        datatype(Float32),
                        dataspace((sizex, sizey, sizez));
                        chunk = (cx, cy, cz),
                    )
                end
            end
        end

        # Create datasets for the prognostic variables.
        if prepare_restart || :rhop in output_variables
            create_dataset(
                file,
                "rhop",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if :u in output_variables
            create_dataset(
                file,
                "u",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if prepare_restart || :us in output_variables
            create_dataset(
                file,
                "us",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if :v in output_variables
            create_dataset(
                file,
                "v",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if prepare_restart || :vs in output_variables
            create_dataset(
                file,
                "vs",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if :w in output_variables
            create_dataset(
                file,
                "w",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if :ws in output_variables
            create_dataset(
                file,
                "ws",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if :wtfc in output_variables
            create_dataset(
                file,
                "wtfc",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if prepare_restart || :wstfc in output_variables
            create_dataset(
                file,
                "wstfc",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if :thetap in output_variables
            create_dataset(
                file,
                "thetap",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end
        if prepare_restart || :pip in output_variables
            create_dataset(
                file,
                "pip",
                datatype(Float32),
                dataspace((sizex, sizey, sizez, 0), (sizex, sizey, sizez, -1));
                chunk = (cx, cy, cz, ct),
            )
        end

        if !(typeof(state.namelists.tracer.tracersetup) <: NoTracer)
            for field in fieldnames(TracerPredictands)
                create_dataset(
                    file,
                    string(field),
                    datatype(Float32),
                    dataspace(
                        (sizex, sizey, sizez, 0),
                        (sizex, sizey, sizez, -1),
                    );
                    chunk = (cx, cy, cz, ct),
                )
            end
        end

        if !(typeof(state.namelists.ice.icesetup) <: NoIce)
            for field in fieldnames(IcePredictands)
                create_dataset(
                    file,
                    string(field),
                    datatype(Float32),
                    dataspace(
                        (sizex, sizey, sizez, 0),
                        (sizex, sizey, sizez, -1),
                    );
                    chunk = (cx, cy, cz, ct),
                )
            end
        end

        if !(typeof(state.namelists.turbulence.turbulencesetup) <: NoTurbulence)
            for field in fieldnames(TurbulencePredictands)
                create_dataset(
                    file,
                    string(field),
                    datatype(Float32),
                    dataspace(
                        (sizex, sizey, sizez, 0),
                        (sizex, sizey, sizez, -1),
                    );
                    chunk = (cx, cy, cz, ct),
                )
            end
        end

        # Create datasets for WKB variables.
        if typeof(testcase) <: AbstractWKBTestCase

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
                            (nray_max, sizex, sizey, sizez + 1, 0),
                            (nray_max, sizex, sizey, sizez + 1, -1),
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
                            (sizex, sizey, sizez, 0),
                            (sizex, sizey, sizez, -1),
                        );
                        chunk = (cx, cy, cz, ct),
                    )
                end
            end
        end

        # Return.
        return
    end

    # Return.
    return
end
