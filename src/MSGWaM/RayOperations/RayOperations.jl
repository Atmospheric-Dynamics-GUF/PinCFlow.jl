module RayOperations

using StaticArrays
using ...Types

include("MergedRays.jl")

include("check_rays.jl")
include("compute_intrinsic_frequency.jl")
include("compute_merge_index.jl")
include("compute_spectral_bounds.jl")
include("copy_rays!.jl")
include("get_physical_extent.jl")
include("get_physical_position.jl")
include("get_spectral_extent.jl")
include("get_spectral_position.jl")
include("get_sufaces.jl")
include("merge_wave_action.jl")
include("remove_rays!.jl")

export MergedRays,
    check_rays,
    compute_intrinsic_frequency,
    compute_merge_index,
    compute_spectral_bounds,
    copy_rays!,
    get_physical_extent,
    get_physical_position,
    get_spectral_extent,
    get_spectral_position,
    get_sufaces,
    merge_wave_action,
    remove_rays!

end