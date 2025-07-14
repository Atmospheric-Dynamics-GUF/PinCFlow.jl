module MeanFlowEffect

using LinearAlgebra
using ...Types
using ...Boundaries
using ..Interpolation
using ..RayUpdate

struct UCHI <: AbstractVariable end 
struct VCHI <: AbstractVariable end 
struct WCHI <: AbstractVariable end 

include("apply_blocked_layer_scheme!.jl")
include("apply_shapiro_filter!.jl")
include("compute_gw_integrals!.jl")
include("compute_gw_tendencies!.jl")
include("compute_horizontal_cell_indices.jl")
include("compute_mean_flow_effect!.jl")
include("smooth_gw_tendencies!.jl")
include("leading_order_tracer_fluxes.jl")
include("compute_leading_order_tracer_fluxes!.jl")

export compute_mean_flow_effect!

end
