using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("triad_wave_propagation2.h5") do data
        plot_output(
            "examples/results/triad_wave_propagation2.svg",
            data,
            ("wavespectrum", 2, 1, 110, 19);
            time_unit = "min",
        )
        return
    end
end





