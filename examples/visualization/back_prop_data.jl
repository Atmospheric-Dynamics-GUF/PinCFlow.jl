using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

function mean_nonzero(A::AbstractArray)::Float64
    s = zero(eltype(A))
    n = 0
    @inbounds for x in A
        if x != 0
            s += x
            n += 1
        end
    end
    return n > 0 ? s / n : zero(s)
end

time_slice = 7
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("back_prop_1.h5") do data
        plot_output(
            "examples/results/back_prop_1.svg",
            data,
            ("zr", 2, 1, 90, time_slice);
            time_unit = "min",
        )
        # position of the wave mode at the initial time
        z_pos = data["zr"][:, 3, 1, :, time_slice]
        wave_number = data["mr"][:, 3, 1, :, time_slice]
        mean_pos = mean_nonzero(z_pos)
        mean_wave_number = mean_nonzero(wave_number)
        println((mean_pos, mean_wave_number))
        return
    end
end
