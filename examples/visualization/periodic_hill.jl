# examples/visualization/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("periodic_hill.h5") do data
    plot_contours(
        "examples/results/periodic_hill.svg",
        data,
        ("w", 1, 1, 1, 2);
    )
    return
end
