module PinCFlowMakieExt

using HDF5
using CairoMakie
using LaTeXStrings
using PinCFlow: @ivy

import PinCFlow: plot_output, set_visualization_theme!, symmetric_contours

include("add_contour_plot!.jl")
include("add_scatter_plot!.jl")

include("plot_output.jl")
include("set_visualization_theme!.jl")
include("symmetric_contours.jl")

end
