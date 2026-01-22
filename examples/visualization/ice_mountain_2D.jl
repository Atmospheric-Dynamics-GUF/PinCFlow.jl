# examples/visualization/ice_mountain_2D.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow


h5open("../ice_mountain_wave.h5") do data
    plot_contours(
        "examples/results/ice_mountain_wave.svg",
        data,
        "q",
        (20, 1, 40, 2);
    )
    return
end
