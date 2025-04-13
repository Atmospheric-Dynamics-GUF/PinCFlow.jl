module BoundaryRays

using ...Types
using ...Boundaries
using ..RayOperations

include("set_boundary_rays!.jl")
include("set_meridional_boundary_rays!.jl")
include("set_meridional_halo_rays!.jl")
include("set_vertical_boundary_rays!.jl")
include("set_zonal_boundary_rays!.jl")
include("set_zonal_halo_rays!.jl")

export set_boundary_rays!,
    set_meridional_boundary_rays!,
    set_vertical_boundary_rays!,
    set_zonal_boundary_rays!

end
