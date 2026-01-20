# examples/visualization/wkb_mountain_wave.jl

using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

include("/home/b/b381734/PF/pinc/ext/PinCFlowMakieExt/plot_output.jl")

h5open("/home/b/b381734/spr/tjl23/ice_mountain_wave.h5") do data
    plot_output(
        "./exp/visualization/test.svg",
        data,
        ("nr", 10, 1, 10, 5);
        time_unit = "s"
    )
    return
end
