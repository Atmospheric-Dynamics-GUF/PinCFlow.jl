using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wkb_wave_propagation_new.h5") do data
        plot_output(
            "examples/results/wkb_wave_propagation_new.svg",
            data,
            ("wavespectrum", 2, 1, 16, 10);
            time_unit = "min",
        )
        return
    end
end





