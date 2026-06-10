"""
```julia
ensemble(simulation::Function, parameters::NamedTuple, keywords::NamedTuple)
```

Run `simulation` in an ensemble.

# Arguments

  - `simulation`: Function to be run in the ensemble. For each key in `parameters` and `keywords`, `simulation` must have a matching keyword argument. In addition, it must have the keyword argument `base_comm`.

  - `parameters`: Keyword arguments for `simulation`, which have different values for different ensemble members. Each entry of `parameters` must be a tuple of ensemble values for the keyword argument represented by the key. One of the keys must be `:output_file`.

  - `keywords`: Keyword arguments for `simulation`, which have the same values for all ensemble members.

  - `delay`: Delay (in seconds) before the first exception is rethrown.
"""
function ensemble end

function ensemble(
    simulation::Function,
    parameters::NamedTuple,
    keywords::NamedTuple;
    delay::Real = 0,
)

    # Check if the output file is one of the parameters.
    !(:output_file in keys(parameters)) &&
        error("The parameter output_file must be assigned in ensembles!")

    (; output_file) = parameters

    # Check if there are duplicate output files.
    for entry in output_file
        sum(1 for other_entry in output_file if other_entry == entry) != 1 &&
            error("There are duplicate output files!")
    end

    # Initialize MPI and get the rank and split color.
    MPI.Init()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    color = rank % length(parameters[keys(parameters)[1]]) + 1
    base_comm = MPI.Comm_split(MPI.COMM_WORLD, color, rank)
    child_rank = MPI.Comm_rank(base_comm)

    # Run the simulations.
    reduce_exceptions(
        MPI.COMM_WORLD;
        delay,
        info = "Ensemble member $(color) has thrown the following exception:",
    ) do
        child_rank == 0 && mkpath(dirname(output_file[color]))

        MPI.Barrier(base_comm)

        open(replace(output_file[color], r"\.h5$" => ".log"), "w") do io
            redirect_stdio(; stderr = io, stdout = io) do
                simulation(;
                    NamedTuple(
                        key => parameters[key][color] for
                        key in keys(parameters)
                    )...,
                    keywords...,
                    base_comm,
                )
                return
            end
            return
        end
        return
    end

    return
end
