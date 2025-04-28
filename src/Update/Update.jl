module Update

using MPI
using ..Types

struct Cartesian end
struct TFC end
struct LHS end
struct RHS end

include("apply_unified_sponge!.jl")
include("compute_sponge!.jl")
include("compute_stress_tensor.jl")
include("compute_time_step.jl")
include("compute_vertical_wind.jl")
include("compute_volume_force.jl")
include("modify_compressible_wind!.jl")
include("synchronize_density_fluctuations!.jl")
include("transform.jl")
include("update_buoyancy_frequency!.jl")
include("update!.jl")

export LHS, RHS

export apply_unified_sponge!,
    compute_sponge!,
    compute_stress_tensor,
    compute_time_step,
    compute_vertical_wind,
    modify_compressible_wind!,
    synchronize_density_fluctuations!,
    update_buoyancy_frequency!,
    update!

end
