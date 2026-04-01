"""
```julia
ensemble(
    code::AbstractString,
    ensemble_size::Integer,
    parameters::Vararg{Symbol};
)::Function
```

Create and return an anonymous function, which can be used to run `code` in an ensemble of size `ensemble_size`, with different values for `parameters` between members.

The returned function takes tuples as arguments which contain the values for each of the parameters. These tuples must be as many and in the same order as the parameters, and each tuple must have the length `ensemble_size`.

# Arguments

  - `code`: Code to be run in the ensemble.

  - `ensemble_size`: Number of ensemble members

  - `parameters`: Variables defined in `code`, which are to be assigned different values for each ensemble member.
"""
function ensemble end

function ensemble(
    code::AbstractString,
    ensemble_size::Integer,
    parameters::Vararg{Symbol},
)::Function

    # Check if the output file is one of the parameters.
    !(:output_file in parameters) &&
        error("The parameter output_file must be assigned in ensembles!")

    # Initialize MPI and get the rank and split color.
    MPI.Init()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    color = rank % ensemble_size + 1

    # Modify the variable assignments.
    modified_code = replace_assignments(
        code,
        [name => Symbol("$(name)[$(color)]") for name in parameters]...,
        :base_comm => :base_comm;
        allow_missing_assignments = false,
    )

    # Create a block expression from the modified code.
    expression = Meta.parseall(modified_code)
    expression.head = :block

    # Wrap the block expression in a function expression.
    @gensym failure exception
    expression = quote
        $(Expr(:tuple, parameters...)) -> begin
            $(failure) = false
            try
                open(
                    replace(output_file[$(color)], r"\.h5$" => ".log"),
                    "w",
                ) do file
                    redirect_stdout(file) do
                        base_comm = MPI.Comm_split(
                            MPI.COMM_WORLD,
                            $(color),
                            $(rank),
                        )
                        $(expression)
                        return
                    end
                    return
                end
            catch $(exception)
                $(failure) = true
                println(
                    "Ensemble member $(color) failed with the following exception:",
                )
                println($(exception))
            end

            MPI.Barrier(MPI.COMM_WORLD)
            $(failure) &&
                error("At least one ensemble member had an error!")

            return
        end
    end

    return Core.eval(@__MODULE__, expression)
end
