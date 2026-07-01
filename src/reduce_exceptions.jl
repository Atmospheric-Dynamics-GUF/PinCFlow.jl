"""
```julia
reduce_exceptions(
    operation::Function,
    comm::MPI.Comm;
    info::AbstractString = "",
)
```

Execute `operation`, catch exceptions in it, and either rethrow the first one (if `comm` contains only one process) or print its error message and exit via `MPI.Abort` (if `comm` contains multiple processes).

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
    window = MPI.Win_create([true], comm)
    reporter = [false]
    rank = 0

    try
        operation()
    catch exception
        MPI.Win_lock(window; rank, type = :exclusive)
        MPI.Get!(reporter, window; rank)
        MPI.Win_flush(window; rank)
        MPI.Put!([false], window; rank)
        MPI.Win_unlock(window; rank)

        if reporter[1]
            sleep(delay)
            flush(stderr)
            println(stderr, info)
            if MPI.Comm_size(comm) == 1
                rethrow(exception)
            else
                Base.display_error(stderr, exception, catch_backtrace())
                flush(stderr)
                MPI.Abort(comm, MPI.Comm_rank(comm))
            end
        end
    end

    MPI.Barrier(comm)
    MPI.free(window)

    return
end
