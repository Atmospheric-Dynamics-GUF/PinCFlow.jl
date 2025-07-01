"""
```julia
MPIOperations
```

Module for operations that require communication across MPI processes.

Provides halo exchange functions for maintaining field continuity across
process boundaries and global reduction operations for distributed arrays.

# Key Functions

  - [`set_zonal_halos_of_field!`](@ref): Halo exchange in ``x``-direction.
  - [`set_meridional_halos_of_field!`](@ref): Halo exchange in ``y``-direction.
  - [`set_vertical_halos_of_field!`](@ref): Halo exchange in ``z``-direction.
  - [`compute_global_dot_product`](@ref): Distributed dot product computation.

# Array Support

  - 2D, 3D, and 5D field arrays
  - Custom halo-layer specifications
  - Solid-wall boundary handling in the vertical.
"""
module MPIOperations

using MPI
using LinearAlgebra
using ..Types

include("compute_global_dot_product.jl")
include("set_meridional_halos_of_field!.jl")
include("set_vertical_halos_of_field!.jl")
include("set_zonal_halos_of_field!.jl")

export compute_global_dot_product,
    set_meridional_halos_of_field!,
    set_vertical_halos_of_field!,
    set_zonal_halos_of_field!

end
