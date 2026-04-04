include("../PinCFlowExamples.jl/src/WavePacketTools/WavePacketTools.jl")
include("TestTools/TestTools.jl")

using Test
using PinCFlow
using .WavePacketTools
using .TestTools

const examples_directory = "../PinCFlowExamples.jl/src/"
const update_references = false

for file in readdir(examples_directory)
    test_file = "test_" * file
    if isfile(test_file)
        include(examples_directory * file)
        include(test_file)
    end
end

@testset verbose = true "PinCFlow tests" begin
    test_cold_bubble()
    test_hot_bubble()
    test_mountain_wave()
    test_periodic_hill()
    test_vortex()
    test_wave_packet()
    test_wkb_mountain_wave()
    test_wkb_wave_packet()
end
