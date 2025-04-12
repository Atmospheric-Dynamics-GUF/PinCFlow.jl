module MeanFlowEffect

using ...Types
using ...Boundaries
using ..Interpolation

include("compute_gw_forcing!.jl")
include("compute_gw_heating!.jl")
include("compute_gw_tendencies!.jl")
include("compute_integrals!.jl")
include("compute_mean_flow_effect!.jl")
include("compute_horizontal_cell_indices.jl")
include("set_integrals_to_zero!.jl")
include("smooth_tendencies!.jl")

export calc_meanflow_effect!

end