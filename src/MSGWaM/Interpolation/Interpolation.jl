module Interpolation

using ...Types
using ...Update

struct N2 <: AbstractVariable end
struct DN2DZ <: AbstractVariable end
struct DUDX <: AbstractVariable end
struct DUDY <: AbstractVariable end
struct DUDZ <: AbstractVariable end
struct DVDX <: AbstractVariable end
struct DVDY <: AbstractVariable end
struct DVDZ <: AbstractVariable end

include("get_next_half_level.jl")
include("get_next_level.jl")
include("interpolate_mean_flow.jl")
include("interpolate_sponge.jl")
include("interpolate_stratification.jl")
include("interpolate.jl")

export get_next_half_level,
    get_next_level,
    interpolate_mean_flow,
    interpolate_sponge,
    interpolate_stratification,
    N2

end
