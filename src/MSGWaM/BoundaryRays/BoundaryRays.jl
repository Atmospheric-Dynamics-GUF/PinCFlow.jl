"""
```julia
BoundaryRays
```

Module for enforcing boundary conditions for ray volumes.

Provides functions for configurations that are serial or parallel in any dimension of physical space. Assumes periodicity in the horizontal and solid-wall boundaries in the vertical.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.MPIOperations`](@ref)

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations`](@ref)
"""
module BoundaryRays

using MPI
using ..RayOperations
using ...Types
using ...MPIOperations
using ...Boundaries
using ...PinCFlow

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
