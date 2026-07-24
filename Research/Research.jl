module Research

using PrecompileTools
using Revise
using PinCFlow
using MPI
using HDF5

include("wkb_wp.jl")
include("wp.jl")

@setup_workload begin
    redirect_stdio(; stderr = devnull, stdout = devnull) do
        mktempdir() do directory
            x_size = 3
            y_size = 3
            z_size = 5

            keywords = (
                output_file = directory * "/pincflow_output.h5",
            )

            @compile_workload begin
                wp(; x_size, y_size, z_size, keywords...)
                wkb_wp(; x_size, z_size, keywords...)
            end
            return
        end
        return
    end
end

export wp_1d, wkb_wp_1d, wkb_wp, wp

end
