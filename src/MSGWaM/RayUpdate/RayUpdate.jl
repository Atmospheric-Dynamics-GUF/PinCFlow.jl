module RayUpdate

using MPI
using Statistics
using ...Types
using ..BoundaryRays
using ..Interpolation
using ..RayOperations
using ..RaySources

struct X end
struct Y end
struct Z end
struct XZ end
struct YZ end
struct XYZ end

include("apply_saturation_scheme!.jl")
include("initialize_rays!.jl")
include("merge_rays!.jl")
include("propagate_rays!.jl")
include("shift_rays!.jl")
include("split_rays!.jl")

export apply_saturation_scheme!,
    initialize_rays!, merge_rays!, propagate_rays!, shift_rays!, split_rays!, X, Y, Z, XZ, YZ, XYZ

end
