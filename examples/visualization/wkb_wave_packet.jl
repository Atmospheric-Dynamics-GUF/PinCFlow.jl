using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("multiple_wkb_wave_packet.h5") do data
        plot_output(
            "examples/results/multiple_wkb_wave_packet1.svg",
            data,
            ("wavespectrum", 15, 15, 16, 4);
            time_unit = "s",
            space_unit = "m"
        )
        return
    end
end