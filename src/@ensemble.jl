"""
```julia
@ensemble(code::Expr, assignments::Expr)
```

Run `code` in an ensemble, where each ensemble member uses different assignment replacements.

# Arguments

  - `code`: Expression that evaluates to a string of Julia code.

  - `assignments`: Expression that evaluates to an array or tuple of pairs, where each pair consists of a symbol (representing the variable name) and an array or tuple that contains the corresponding values that are to be used by the ensemble members. Note that for every replacement given in this argument, a corresponding assignment must be present in `code`. Replacements of the variable output_file are required.
"""
macro ensemble end

macro ensemble(code::Expr, assignments::Expr)
    return quote
        local index = findlast(
            [name for (name, value) in $(esc(assignments))] .== :output_file,
        )
        index === nothing &&
            error("The parameter output_file must be assigned in ensembles!")
        local output_file = $(esc(assignments))[index][2]
        for entry in output_file
            sum(1 for other_entry in output_file if other_entry == entry) !=
            1 && error("There are duplicate output files!")
        end

        MPI.Init()
        local rank = MPI.Comm_rank(MPI.COMM_WORLD)
        local ensemble_size = length($(esc(assignments))[1][2])
        local color = rank % ensemble_size + 1

        local modified_code = replace_assignments(
            $(esc(code)),
            [name => value[color] for (name, value) in $(esc(assignments))]...,
            :base_comm => MPI.Comm_split(MPI.COMM_WORLD, color, rank);
            allow_missing_assignments = false,
        )

        local failure = false
        try
            open(replace(output_file[color], r"\.h5$" => ".log"), "w") do file
                redirect_stdout(file) do
                    Core.eval(@__MODULE__, Meta.parseall(modified_code))
                    return
                end
            end
        catch exception
            local failure = true
            println(
                "Ensemble member $(color) failed with the following exception:",
            )
            println(exception)
            println("")
        end

        MPI.Barrier(MPI.COMM_WORLD)

        failure && error("At least one ensemble member had an error!")
    end
end
