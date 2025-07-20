"""
```julia
OutputNamelist{
    A <: Tuple{Vararg{Symbol, <:Integer}},
    B <: Bool,
    C <: Integer,
    D <: AbstractFloat,
    E <: AbstractString,
}
```

Namelist for I/O parameters.

# Fields

  - `output_variables::A`: A tuple of symbols representing the variables that should be written to the output file.
  - `save_ray_volumes::B`: A boolean indicating whether to write ray-volume data.
  - `prepare_restart::B`: A boolean indicating whether to write all variables needed for restart simulations.
  - `restart::B`: A boolean indicating whether to initialize with data from a previous state (as written in `input_file`).
  - `iin::C`: Temporal index in `input_file` at which to read the data to initialize with in restart simulations.
  - `output_steps::B`: If set to `true`, write output every `noutput` time steps.
  - `noutput::C`: Output interval (in indices) if `output_steps == true`.
  - `maxiter::C`: Maximum number of iterations if `output_steps == true`.
  - `outputtimediff::D`: Output interval (in physical time) if `output_steps == false`.
  - `maxtime::D`: Simulation time if `output_steps == false`.
  - `input_file::E`: File from which to read input data in restart simulations.
  - `output_file::E`: File to which output data is written.
"""
struct OutputNamelist{
    A <: Tuple{Vararg{Symbol, <:Integer}},
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
    noutput::C
    maxiter::C
    outputtimediff::D
    maxtime::D
    input_file::E
    output_file::E
end

"""
```julia
OutputNamelist(;
    output_variables::Tuple{Vararg{Symbol, <:Integer}} = (),
    save_ray_volumes::Bool = false,
    prepare_restart::Bool = false,
    restart::Bool = false,
    iin::Integer = -1,
    output_steps::Bool = false,
    noutput::Integer = 1,
    maxiter::Integer = 1,
    outputtimediff::AbstractFloat = 3.6E+3,
    maxtime::AbstractFloat = 3.6E+3,
    input_file::AbstractString = "./pincflow_input.h5",
    output_file::AbstractString = "./pincflow_output.h5",
)
```

Construct an `OutputNamelist` instance with the given keyword arguments as properties.

# Arguments

  - `output_variables`: A tuple of symbols representing the variables that should be written to the output file.
  - `save_ray_volumes`: A boolean indicating whether to write ray-volume data.
  - `prepare_restart`: A boolean indicating whether to write all variables needed for restart simulations.
  - `restart`: A boolean indicating whether to initialize with data from a previous state (as written in `input_file`).
  - `iin`: Temporal index in `input_file` at which to read the data to initialize with in restart simulations.
  - `output_steps`: If set to `true`, write output every `noutput` time steps.
  - `noutput`: Output interval (in indices) if `output_steps == true`.
  - `maxiter`: Maximum number of iterations if `output_steps == true`.
  - `outputtimediff`: Output interval (in physical time) if `output_steps == false`.
  - `maxtime`: Simulation time if `output_steps == false`.
  - `input_file`: File from which to read input data in restart simulations.
  - `output_file`: File to which output data is written.

# Returns

  - `::OutputNamelist`: `OutputNamelist` instance.
"""
function OutputNamelist(;
    output_variables::Tuple{Vararg{Symbol, <:Integer}} = (),
    save_ray_volumes::Bool = false,
    prepare_restart::Bool = false,
    restart::Bool = false,
    iin::Integer = -1,
    output_steps::Bool = false,
    noutput::Integer = 1,
    maxiter::Integer = 1,
    outputtimediff::AbstractFloat = 3.6E+3,
    maxtime::AbstractFloat = 3.6E+3,
    input_file::AbstractString = "./pincflow_input.h5",
    output_file::AbstractString = "./pincflow_output.h5",
)
    return OutputNamelist(
        output_variables,
        save_ray_volumes,
        prepare_restart,
        restart,
        iin,
        output_steps,
        noutput,
        maxiter,
        outputtimediff,
        maxtime,
        input_file,
        output_file,
    )
end
