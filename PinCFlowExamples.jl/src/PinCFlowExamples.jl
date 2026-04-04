module PinCFlowExamples

using MPI
using HDF5
using CairoMakie
using PrecompileTools
using Revise
using PinCFlow

include("WavePacketTools/WavePacketTools.jl")

using .WavePacketTools

include("cold_bubble.jl")
include("hot_bubble.jl")
include("mountain_wave.jl")
include("periodic_hill.jl")
include("vortex.jl")
include("wave_packet.jl")
include("wkb_mountain_wave.jl")
include("wkb_wave_packet.jl")

@setup_workload begin
    @compile_workload begin
        redirect_stdout(devnull) do
            mktempdir() do directory
                x_size = 5
                y_size = 5
                z_size = 5

                npx = 1
                npy = 1
                npz = 1

                output_file = directory * "/pincflow_output.h5"
                output_steps = true

                visualize = false

                cold_bubble(;
                    x_size,
                    z_size,
                    npx,
                    npz,
                    output_file,
                    output_steps,
                    visualize,
                )

                hot_bubble(;
                    x_size,
                    z_size,
                    npx,
                    npz,
                    output_file,
                    output_steps,
                    visualize,
                )

                mountain_wave(;
                    x_size,
                    y_size,
                    z_size,
                    npx,
                    npy,
                    npz,
                    output_file,
                    output_steps,
                    visualize,
                )

                periodic_hill(;
                    x_size,
                    z_size,
                    npx,
                    npz,
                    output_file,
                    output_steps,
                    visualize,
                )

                vortex(;
                    x_size,
                    y_size,
                    npx,
                    npy,
                    output_file,
                    output_steps,
                    visualize,
                )

                wave_packet(;
                    x_size,
                    y_size,
                    z_size,
                    npx,
                    npy,
                    npz,
                    output_file,
                    output_steps,
                    visualize,
                )

                wkb_mountain_wave(;
                    x_size,
                    y_size,
                    z_size,
                    npx,
                    npy,
                    npz,
                    output_file,
                    output_steps,
                    visualize,
                )

                wkb_wave_packet(;
                    x_size,
                    y_size,
                    z_size,
                    npx,
                    npy,
                    npz,
                    output_file,
                    output_steps,
                    visualize,
                )

                return
            end
            return
        end
    end
end

export cold_bubble,
    hot_bubble,
    mountain_wave,
    periodic_hill,
    vortex,
    wave_packet,
    wkb_mountain_wave,
    wkb_wave_packet

end
