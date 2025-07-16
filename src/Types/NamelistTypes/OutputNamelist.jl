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

Namelist for I/O (see constructor for parameter descriptions).
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
    output_variables = (),
    save_ray_volumes = false,
    prepare_restart = false,
    restart = false,
    iin = -1,
    output_steps = false,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 3.6E+3,
    maxtime = 3.6E+3,
    input_file = "./pincflow_input.h5",
    output_file = "./pincflow_output.h5",
)
```

Construct a new OutputNamelist instance, which holds output and input related parameters.

# Arguments

  - `output_variables`: A tuple of symbols representing the output variables which should be written to the output file.
  - `save_ray_volumes`: A boolean indicating whether to write out ray volumes.
  - `prepare_restart`: A boolean indicating whether to prepare an output file for restart.
  - `restart`: A boolean indicating whether to restart from a previous state. If true, restart from a previous state saved in `input_file`.
  - `iin`: An integer representing the time index used for restart.
  - `output_steps`: If true, write output every `noutput` steps.
  - `noutput`: If `output_steps` is true, an integer representing the number of steps between outputs.
  - `maxiter`: An integer representing the maximum number of iterations. Only used if `output_steps` is true.
  - `outputtimediff`: a floating-point number determining the time difference between outputs in seconds. Only used if `output_steps` is false.
  - `maxtime`: A floating-point number representing the maximum simulaton time. Only used if `output_steps` is false.
  - `fancy_namelists`: Not used for now.
  - `input_file`: A string holding the input HDF5 file path used for restart.
  - `output_file`: A string holding the output HDF5 file path.

# Returns

  - `::OutputNamelist`: `OutputNamelist` instance.
"""
function OutputNamelist(;
    output_variables = (),
    save_ray_volumes = false,
    prepare_restart = false,
    restart = false,
    iin = -1,
    output_steps = false,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 3.6E+3,
    maxtime = 3.6E+3,
    input_file = "./pincflow_input.h5",
    output_file = "./pincflow_output.h5",
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
