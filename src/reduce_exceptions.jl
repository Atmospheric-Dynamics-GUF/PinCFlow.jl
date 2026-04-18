"""
```julia
reduce_exceptions(
    operation::Function,
    comm::MPI.Comm;
    info::AbstractString = "",
)
```

Execute `operation`, catch exceptions in it, and rethrow the exception in the fastest process.

# Arguments

  - `operation`: Function which takes no arguments and returns `nothing`.

  - `comm`: MPI communicator containing all processes participating in `operation`.

  - `delay`: Delay (in seconds) before the first exception is rethrown.

  - `info`: String to print just before the caught exception is rethrown.
"""
function reduce_exceptions end

function reduce_exceptions(
    operation::Function,
    comm::MPI.Comm;
    delay::Real = 0,
    info::AbstractString = "",
)
    lock = MPI.bcast(tempname(), comm)

    try
        operation()
    catch exception
        reporter = false
        try
            mkdir(lock)
            reporter = true
        catch
        end

        if reporter
            sleep(delay)
            flush(stdout)
            println(info)
            rethrow(exception)
        end
    end

    MPI.Barrier(comm)

    return
end
