"""
```julia
OutputNamelist{
    A <: Tuple{Vararg{Symbol}},
    B <: Bool,
    C <: Integer,
    D <: AbstractFloat,
    E <: AbstractString,
}
```

Namelist for I/O parameters.

```julia
OutputNamelist(;
    output_variables::Tuple{Vararg{Symbol}} = (),
    save_ray_volumes::Bool = false,
    prepare_restart::Bool = false,
    restart::Bool = false,
    iin::Integer = -1,
    output_steps::Bool = false,
    nout::Integer = 1,
    iterations::Integer = 1,
    output_interval::AbstractFloat = 3.6E+3,
    tmax::AbstractFloat = 3.6E+3,
    input_file::AbstractString = "./pincflow_input.h5",
    output_file::AbstractString = "./pincflow_output.h5",
)::OutputNamelist
```

Construct an `OutputNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `output_variables::A`: A tuple of symbols representing the variables that should be written to the output file.

  - `save_ray_volumes::B`: A boolean indicating whether to write ray-volume data.

  - `prepare_restart::B`: A boolean indicating whether to write all variables needed for restart simulations.

  - `restart::B`: A boolean indicating whether to initialize with data from a previous state (as written in `input_file`).

  - `iin::C`: Temporal index in `input_file` at which to read the data to initialize with in restart simulations.

  - `output_steps::B`: If set to `true`, write output every `nout` time steps.

  - `nout::C`: Output interval (in indices) if `output_steps == true`.

  - `iterations::C`: Maximum number of iterations if `output_steps == true`.

  - `output_interval::D`: Output interval (in physical time) if `output_steps == false`.

  - `tmax::D`: Simulation time if `output_steps == false`.

  - `input_file::E`: File from which to read input data in restart simulations.

  - `output_file::E`: File to which output data is written.
"""
struct OutputNamelist{
    A <: Tuple{Vararg{Symbol}},
    B <: Bool,
    C <: Integer,
    D <: AbstractFloat,
    E <: AbstractString,
}
    output_variables::A
    save_ray_volumes::B
    prepare_restart::B
    restart::B
    iin::C
    output_steps::B
    nout::C
    iterations::C
    output_interval::D
    tmax::D
    input_file::E
    output_file::E
end

function OutputNamelist(;
    output_variables::Tuple{Vararg{Symbol}} = (),
    save_ray_volumes::Bool = false,
    prepare_restart::Bool = false,
    restart::Bool = false,
    iin::Integer = -1,
    output_steps::Bool = false,
    nout::Integer = 1,
    iterations::Integer = 1,
    output_interval::AbstractFloat = 3.6E+3,
    tmax::AbstractFloat = 3.6E+3,
    input_file::AbstractString = "./pincflow_input.h5",
    output_file::AbstractString = "./pincflow_output.h5",
)::OutputNamelist
    return OutputNamelist(
        output_variables,
        save_ray_volumes,
        prepare_restart,
        restart,
        iin,
        output_steps,
        nout,
        iterations,
        output_interval,
        tmax,
        input_file,
        output_file,
    )
end
