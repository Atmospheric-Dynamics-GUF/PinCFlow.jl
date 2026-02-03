using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow




#if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("triad_wave_propagation_250.h5") do data
        my_plot_output(
            "examples/results/triad_wave_propagation_250_time14.svg",
            data,
            ("wavespectrum", 1, 1, 60, 65);
            space_unit = "m",
            time_unit = "min",
        )
        return
    end
#end





