using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("resolved_wave_prop_without_shear.h5") do data
        plot_output(
            "examples/results/resolved_wave_prop_without_shear.svg",
            data,
            ("w", 33, 1, 160, 1);
            time_unit = "min",
        )
        return
    end
end
