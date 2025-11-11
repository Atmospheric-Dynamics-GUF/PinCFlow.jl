# examples/visualization/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("prop_perturb.h5") do data
    plot_contours(
       "examples/results/prop_perturb.svg",
        data,
        "w",
        (10, 10, 1, 1);
        label = L"w\,[\mathrm{m\,s^{-1}}]",
    )
    return
end





