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
    fancy_namelists::B
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
    fancy_namelists = true,
    input_file = "./pincflow_input.h5",
    output_file = "./pincflow_output.h5",
)
```

Construct a new OutputNamelist instance, which holds output and input related parameters.

# Arguments

  - `output_variables::A`: A tuple of symbols representing the output variables which should be written to the output file.
  - `save_ray_volumes::B`: A boolean indicating whether to write out ray volumes.
  - `prepare_restart::B`: A boolean indicating whether to prepare an output file for restart.
  - `restart::B`: A boolean indicating whether to restart from a previous state. If true, restart from a previous state saved in `input_file`.
  - `iin::C`: An integer representing the time index used for restart.
  - `output_steps::B`: If true, write output every `noutput` steps.
  - `noutput::C`: If `output_steps` is true, an integer representing the number of steps between outputs.
  - `maxiter::C`: An integer representing the maximum number of iterations. Only used if `output_steps` is true.
  - `outputtimediff::D`: a floating-point number determining the time difference between outputs in seconds. Only used if `output_steps` is false.
  - `maxtime::D`: A floating-point number representing the maximum simulaton time. Only used if `output_steps` is false.
  - `fancy_namelists::B`: Not used for now.
  - `input_file::E`: A string holding the input HDF5 file path used for restart.
  - `output_file::E`: A string holding the output HDF5 file path.
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
    fancy_namelists = true,
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
        fancy_namelists,
        input_file,
        output_file,
    )
end
