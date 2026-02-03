# examples/visualization/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("resolved_flow_test_1.h5") do data
    plot_contours(
       "examples/results/resolved_flow_test_1.svg",
        data,
        "w",
        (2, 2, 2, 4);
        label = L"w\,[\mathrm{m\,s^{-1}}]",
    )
    return
end





