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
