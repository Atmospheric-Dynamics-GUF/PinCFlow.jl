# Package extension for adding PythonPlot-based features to PinCFlow.jl.
module PinCFlowPythonPlotExt

using PythonPlot

using PinCFlow: @ivy

include("PythonPlotExt/set_plot_style.jl")
include("PythonPlotExt/symmetric_contours.jl")

end