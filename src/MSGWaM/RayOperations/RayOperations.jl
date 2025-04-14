module RayOperations

using StaticArrays
using ...Types
using ..Interpolation

include("MergedRays.jl")

include("check_rays.jl")
include("compute_intrinsic_frequency.jl")
include("compute_merge_index.jl")
include("compute_saturation_integrals.jl")
include("compute_spectral_bounds.jl")
include("copy_rays!.jl")
include("get_physical_extent.jl")
include("get_physical_position.jl")
include("get_spectral_extent.jl")
include("get_spectral_position.jl")
include("get_surfaces.jl")
include("merge_wave_action.jl")
include("remove_rays!.jl")

export MergedRays,
    check_rays,
    compute_intrinsic_frequency,
    compute_merge_index,
    compute_saturation_integrals,
    compute_spectral_bounds,
    copy_rays!,
    get_physical_extent,
    get_physical_position,
    get_spectral_extent,
    get_spectral_position,
    get_surfaces,
    merge_wave_action,
    remove_rays!

end
