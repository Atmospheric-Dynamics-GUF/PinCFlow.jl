"""
```julia
OutputNamelist
```

Namelist for I/O parameters.

```julia
OutputNamelist(;
    output_variables::Vector{Symbol} = Symbol[],
    save_ray_volumes::Bool = false,
    prepare_restart::Bool = false,
    restart::Bool = false,
    iin::Integer = -1,
    output_steps::Bool = false,
    nout::Integer = 1,
    iterations::Integer = 1,
    output_interval::Real = 3.6E+3,
    tmax::Real = 3.6E+3,
    input_file::AbstractString = "./pincflow_input.h5",
    output_file::AbstractString = "./pincflow_output.h5",
)::OutputNamelist
```

Construct an `OutputNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `output_variables::Vector{Symbol}`: A vector of symbols representing the variables that should be written to the output file.

  - `save_ray_volumes::Bool`: A boolean indicating whether to write ray-volume data.

  - `prepare_restart::Bool`: A boolean indicating whether to write all variables needed for restart simulations.

  - `restart::Bool`: A boolean indicating whether to initialize with data from a previous state (as written in `input_file`).

  - `iin::Int`: Temporal index in `input_file` at which to read the data to initialize with in restart simulations. If it's set to the default value `-1`, the data will be read at the last index.

  - `output_steps::Bool`: If set to `true`, write output every `nout` time steps.

  - `nout::Int`: Output interval (in indices) if `output_steps == true`.

  - `iterations::Int`: Maximum number of iterations if `output_steps == true`.

  - `output_interval::Float64`: Output interval (in physical time) if `output_steps == false`.

  - `tmax::Float64`: Simulation time if `output_steps == false`.

  - `input_file::String`: File from which to read input data in restart simulations.

  - `output_file::String`: File to which output data is written.
"""
struct OutputNamelist
    output_variables::Vector{Symbol}
    save_ray_volumes::Bool
    prepare_restart::Bool
    restart::Bool
    iin::Int
    output_steps::Bool
    nout::Int
    iterations::Int
    output_interval::Float64
    tmax::Float64
    input_file::String
    output_file::String
end

function OutputNamelist(;
    output_variables::Vector{Symbol} = Symbol[],
    save_ray_volumes::Bool = false,
    prepare_restart::Bool = false,
    restart::Bool = false,
    iin::Integer = -1,
    output_steps::Bool = false,
    nout::Integer = 1,
    iterations::Integer = 1,
    output_interval::Real = 3.6E+3,
    tmax::Real = 3.6E+3,
    input_file::AbstractString = "./pincflow_input.h5",
    output_file::AbstractString = "./pincflow_output.h5",
)::OutputNamelist
    return OutputNamelist(
        output_variables,
        save_ray_volumes,
        prepare_restart,
        restart,
        Int(iin),
        output_steps,
        Int(nout),
        Int(iterations),
        Float64(output_interval),
        Float64(tmax),
        string(input_file),
        string(output_file),
    )
end
