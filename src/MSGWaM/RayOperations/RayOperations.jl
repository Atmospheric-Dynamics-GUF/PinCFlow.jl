"""
```julia
RayOperations
```

Module for various ray-volume operations needed throughout `PinCFlow.MSGWaM`.

# Scopes

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation`](@ref)

  - [`PinCFlow.Types`](@ref)

# Exports

  - [`PinCFlow.MSGWaM.RayOperations.check_rays`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_intrinsic_frequency`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_merge_index`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_saturation_integrals`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_spectral_bounds`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_turbulent_velocity`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_wave_action_integral`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_physical_extent`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_physical_position`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_spectral_extent`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_spectral_position`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_surfaces`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.remove_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.update_merged_rays!`](@ref)
"""
module RayOperations

using ..Interpolation
using ...Macros
using ...Types

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
include("compute_turbulent_velocity.jl")

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
    remove_rays!,
    update_merged_rays!,
    compute_turbulent_velocity

end
