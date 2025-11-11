# examples/visualization/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("wkb_prop_perturb.h5") do data
    plot_contours(
       "examples/results/wkb_prop_perturb.svg",
        data,
        "dudt",
        (2, 2, 2, 2);
        label = L"w\,[\mathrm{m\,s^{-1}}]",
    )
    return
end





