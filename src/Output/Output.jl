"""
```julia
Output
```

Module for I/O of simulation data.

Provides functions for writing the model state and initializing the model with data from a previous simulation.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Update`](@ref)
"""
module Output

using Dates
using MPI
using HDF5
using ..Types
using ..Update
using ..PinCFlow

include("create_output.jl")
include("read_input!.jl")
include("write_output.jl")

export create_output, read_input!, write_output

end
