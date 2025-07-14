module Update

using MPI
using ..Types
using ..Boundaries

struct Cartesian end
struct Transformed end
struct LHS end
struct RHS end
struct IceUpdateAdv end
struct IceUpdatePhy end
struct IceUpdateBoth end

include("apply_unified_sponge!.jl")
include("compute_compressible_buoyancy_factor.jl")
include("compute_compressible_wind_factor.jl")
include("compute_sponge!.jl")
include("compute_stress_tensor.jl")
include("compute_vertical_wind.jl")
include("compute_volume_force.jl")
include("transform.jl")
include("update!.jl")

export LHS, RHS, IceUpdatePhy

export apply_unified_sponge!,
    compute_compressible_buoyancy_factor,
    compute_compressible_wind_factor,
    compute_sponge!,
    compute_stress_tensor,
    compute_vertical_wind,
    compute_volume_force,
    update!

end
