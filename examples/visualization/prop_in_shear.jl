# examples/visualization/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("prop_without_shear.h5") do data
    plot_contours(
       "examples/results/prop_without_shear.svg",
        data,
        "w",
        (2, 2, 2, 4);
        label = L"w\,[\mathrm{m\,s^{-1}}]",
    )
    return
end





