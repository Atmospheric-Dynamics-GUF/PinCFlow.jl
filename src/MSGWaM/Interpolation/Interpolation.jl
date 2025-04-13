module Interpolation

using ...Types
using ...Update

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
    interpolate_stratification

end
