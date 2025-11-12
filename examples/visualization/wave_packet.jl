# examples/visualization/wave_packet.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("wave_packet.h5") do data
    plot_output(
        "examples/results/wave_packet.svg",
        data,
        ("u", 20, 20, 30, 2),
        ("v", 20, 20, 30, 2),
        ("w", 20, 20, 30, 2);
        time_unit = "min",
    )
    return
end
