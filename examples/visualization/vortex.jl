# examples/visualization/vortex.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("vortex.h5") do data
    plot_contours("examples/results/vortex.svg", data, ("chi", 1, 1, 1, 2);)
    return
end
