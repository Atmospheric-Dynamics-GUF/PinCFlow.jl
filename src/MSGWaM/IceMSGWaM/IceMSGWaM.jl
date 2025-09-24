module IceMSGWaM

using ...Types
using ..Interpolation

include("../MeanFlowEffect/compute_horizontal_cell_indices.jl")
include("compute_msgwam_ice!.jl")

export compute_msgwam_ice!

end