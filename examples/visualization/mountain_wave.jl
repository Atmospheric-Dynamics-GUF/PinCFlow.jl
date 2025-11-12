# examples/visualization/mountain_wave.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("mountain_wave.h5") do data
    plot_output(
        "examples/results/mountain_wave.svg",
        data,
        ("w", 20, 20, 10, 2);
    )
    return
end
