# Package extension for adding PythonPlot-based features to PinCFlow.jl.
module PinCFlowPythonPlotExt

using PythonPlot

using PinCFlow: @ivy

import PinCFlow: set_plot_style, symmetric_contours

include("PythonPlotExt/set_plot_style.jl")
include("PythonPlotExt/symmetric_contours.jl")

end