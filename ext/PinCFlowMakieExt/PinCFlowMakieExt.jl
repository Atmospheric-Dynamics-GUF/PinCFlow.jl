module PinCFlowMakieExt

using HDF5
using CairoMakie
using PinCFlow: @ivy

import PinCFlow: plot_output, set_visualization_theme!, symmetric_contours

include("plot_output.jl")
include("set_visualization_theme!.jl")
include("symmetric_contours.jl")

end
