# examples/visualization/hot_bubble.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("hot_bubble.h5") do data
    plot_output(
        "examples/results/hot_bubble.svg",
        data,
        ("thetap", 1, 1, 1, 2);
    )
    return
end
