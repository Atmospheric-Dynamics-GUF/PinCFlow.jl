"""
```julia
RayOperations
```

Module for various ray-volume operations needed throughout `PinCFlow.MSGWaM`.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation`](@ref)
"""
module RayOperations

using ..Interpolation
using ...Types
using ...PinCFlow

include("check_rays.jl")
include("compute_intrinsic_frequency.jl")
include("compute_merge_index.jl")
include("compute_saturation_integrals.jl")
include("compute_spectral_bounds.jl")
include("compute_wave_action_integral.jl")
include("copy_rays!.jl")
include("get_physical_extent.jl")
include("get_physical_position.jl")
include("get_spectral_extent.jl")
include("get_spectral_position.jl")
include("get_surfaces.jl")
include("remove_rays!.jl")
include("update_merged_rays!.jl")

export check_rays,
    compute_intrinsic_frequency,
    compute_merge_index,
    compute_saturation_integrals,
    compute_spectral_bounds,
    compute_wave_action_integral,
    copy_rays!,
    get_physical_extent,
    get_physical_position,
    get_spectral_extent,
    get_spectral_position,
    get_surfaces,
    merge_wave_action,
    remove_rays!,
    update_merged_rays!

end
