# examples/visualization/cold_bubble.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("cold_bubble.h5") do data
    plot_output(
        "examples/results/cold_bubble.svg",
        data,
        ("thetap", 1, 1, 1, 2);
    )
    return
end
