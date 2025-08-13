"""
```julia
MSGWaM
```

3D transient implementation of MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.BoundaryRays`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect`](@ref)

# External links

 1. [Muraschko et al. (2014)](https://doi.org/10.1002/qj.2381)

 1. [Boeloeni et al. (2016)](https://doi.org/10.1175/JAS-D-16-0069.1)

 1. [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1)

 1. [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1)

 1. [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1)
"""
module MSGWaM

include("Interpolation/Interpolation.jl")
include("RayOperations/RayOperations.jl")
include("RaySources/RaySources.jl")
include("BoundaryRays/BoundaryRays.jl")
include("RayUpdate/RayUpdate.jl")
include("MeanFlowEffect/MeanFlowEffect.jl")

using .BoundaryRays
using .RayUpdate
using .MeanFlowEffect

export apply_saturation_scheme!,
    compute_mean_flow_effect!,
    initialize_rays!,
    merge_rays!,
    propagate_rays!,
    set_boundary_rays!,
    shift_rays!,
    split_rays!

end
