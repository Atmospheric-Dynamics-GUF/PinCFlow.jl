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

x_slice = 1
y_slize = 1
z_slice = 110
time_slice = 51
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("back_prop_pos.h5") do data
        plot_output(
            "examples/results/back_prop_pos.svg",
            data,
            ("nr", 1, 1, 110, time_slice);
            time_unit = "min",
        )
        # position of the wave mode at the initial time
        #to find the average hieght and wave numbers
        z_pos = data["zr"][:, 1, 1, :, time_slice]
        wave_number = data["mr"][:, 1, 1, :, time_slice]
        mean_pos = mean_nonzero(z_pos)
        mean_wave_number = mean_nonzero(wave_number)
        println((mean_pos, mean_wave_number))

        #to find the postion and wavenumber towards the maximum nr
        nr_Z = data["nr"][:, 1, 1, :, time_slice]
        peak_nr = argmax(nr_Z)
        (max_alpha, max_z) = Tuple(peak_nr)
        z_pos = data["zr"][max_alpha, 1, 1, max_z, time_slice]
        wave_number = data["mr"][max_alpha, 1, 1, max_z, time_slice]
        println((mean_pos, mean_wave_number))

        return
    end
end
