"""
```julia
BoundaryRays
```

Module for setting the ray volumes in boundary cells.

Provides functions for configurations that are serial/parallel in any dimension of physical space. Assumes periodicity in the horizontal and solid-wall boundaries in the vertical.
"""
module BoundaryRays

using MPI
using ...Types
using ...MPIOperations
using ...Boundaries
using ..RayOperations

include("set_boundary_rays!.jl")
include("set_meridional_boundary_rays!.jl")
include("set_meridional_halo_rays!.jl")
include("set_vertical_boundary_rays!.jl")
include("set_vertical_halo_rays!.jl")
include("set_zonal_boundary_rays!.jl")
include("set_zonal_halo_rays!.jl")

export set_boundary_rays!,
    set_meridional_boundary_rays!,
    set_vertical_boundary_rays!,
    set_zonal_boundary_rays!

end
