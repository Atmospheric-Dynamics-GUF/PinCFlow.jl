"""
```julia
BoundaryRays
```

Module for enforcing boundary conditions for ray volumes.

Provides functions for configurations that are serial or parallel in any dimension of physical space. Assumes periodicity in the horizontal and solid-wall boundaries in the vertical.

# Scopes

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.MPIOperations`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations`](@ref)

  - [`PinCFlow.Types`](@ref)

# Exports

  - [`PinCFlow.MSGWaM.BoundaryRays.set_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_meridional_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_vertical_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_zonal_boundary_rays!`](@ref)
"""
module BoundaryRays

using MPI
using ..RayOperations
using ...Macros
using ...Types
using ...MPIOperations
using ...Boundaries

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
