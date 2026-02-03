using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

#if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("resolved_flow_test_1.h5") do data
        plot_output(
            "examples/results/resolved_flow_test_1.svg",
            data,
            ("w", 40, 1, 1200, 1);
            time_unit = "min",
        )
        return
    end
#end
