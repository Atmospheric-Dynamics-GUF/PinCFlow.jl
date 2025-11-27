using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

slice = 20

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wave_packet1.h5") do data
        plot_output(
            "examples/results/wave_packet1.svg",
            data,
            ("u", slice, 1, 2*slice, 1),
            ("v", slice, 1, 2*slice, 1),
            ("w", slice, 1, 2*slice, 1);
            time_unit = "min",
        )
        return
    end
end