using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

slice = 20

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wave_packet.h5") do data
        plot_output(
            "examples/results/wave_packet.svg",
            data,
            ("u", slice, slice, 2*slice, 2),
            ("v", slice, slice, 2*slice, 2),
            ("w", slice, slice, 2*slice, 2);
            time_unit = "min",
        )
        return
    end
end