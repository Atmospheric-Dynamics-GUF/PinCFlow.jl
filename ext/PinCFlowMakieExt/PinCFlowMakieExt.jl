module PinCFlowMakieExt

using CairoMakie
using PinCFlow: @ivy

import PinCFlow: set_visualization_theme!, symmetric_contours

include("set_visualization_theme!.jl")
include("symmetric_contours.jl")

end
