"""
```julia
MPIOperations
```

Module for operations that require communication between MPI processes.

Provides halo exchange functions for maintaining field continuity across
process boundaries, as well as global reduction operations.

# See also

  - [`PinCFlow.Types`](@ref)
"""
module MPIOperations

using MPI
using ..Types
using ..PinCFlow

include("compute_global_dot_product.jl")
include("set_meridional_halos_of_field!.jl")
include("set_vertical_halos_of_field!.jl")
include("set_zonal_halos_of_field!.jl")

export compute_global_dot_product,
    set_meridional_halos_of_field!,
    set_vertical_halos_of_field!,
    set_zonal_halos_of_field!

end
