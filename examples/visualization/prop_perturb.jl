using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("triad_wave_propagation.h5") do data
        plot_output(
            "examples/results/triad_wave_propagation.svg",
            data,
            ("wavespectrum", 2, 1, 8, 11);
            time_unit = "min",
        )
        return
    end
end





