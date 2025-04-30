module Update

using MPI
using ..Types

struct Cartesian end
struct Transformed end
struct LHS end
struct RHS end

include("apply_unified_sponge!.jl")
include("compute_sponge!.jl")
include("compute_stress_tensor.jl")
include("compute_vertical_wind.jl")
include("compute_volume_force.jl")
include("transform.jl")
include("update!.jl")

export LHS, RHS

export apply_unified_sponge!,
    compute_sponge!,
    compute_stress_tensor,
    compute_vertical_wind,
    compute_volume_force,
    update!

end
