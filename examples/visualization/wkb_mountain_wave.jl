# examples/visualization/wkb_mountain_wave.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("wkb_mountain_wave.h5") do data
    plot_output(
        "examples/results/wkb_mountain_wave.svg",
        data,
        ("ar", 20, 20, 10, 2);
    )
    return
end
