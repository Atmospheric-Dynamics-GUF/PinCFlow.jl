module MSGWaM

include("Interpolation/Interpolation.jl")
include("RayOperations/RayOperations.jl")
include("RaySources/RaySources.jl")
include("BoundaryRays/BoundaryRays.jl")
include("RayUpdate/RayUpdate.jl")
include("MeanFlowEffect/MeanFlowEffect.jl")

using .RayUpdate
using .MeanFlowEffect

export apply_saturation_scheme!,
    initialize_rays!, merge_rays!, propagate_rays!, shift_rays!, split_rays!

end