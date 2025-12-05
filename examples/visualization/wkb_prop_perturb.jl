# examples/visualization/periodic_hill.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("wkb_wave_packet.h5") do data
    plot_contours(
       "examples/results/wkb_wave_packet.svg",
        data,
        "w",
        (5, 5, 5, 4);
        label = L"w\,[\mathrm{m\,s^{-1}}]",
    )
    return
end





