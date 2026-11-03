include("TestTools/TestTools.jl")

using Aqua
using JuliaFormatter
using MacroTools
using PinCFlow
using Test
using .TestTools

using MacroTools: postwalk

const update = "--update" in ARGS
const keywords = (prepare_restart = true, visualize = false)

include("test_cold_bubble.jl")
include("test_hot_bubble.jl")
include("test_mountain_wave.jl")
include("test_periodic_hill.jl")
include("test_vortex.jl")
include("test_wave_packet.jl")
include("test_wkb_mountain_wave.jl")
include("test_wkb_wave_packet.jl")

include("test_formatting.jl")

include("test_files.jl")

include("test_ascii.jl")
include("test_julia_validity.jl")
include("test_literal_checks.jl")
include("test_long_named_function_definitions.jl")
include("test_type_checks.jl")

@testset verbose = true "PinCFlow tests" begin
    @testset verbose = true "Examples" begin
        test_cold_bubble()
        test_hot_bubble()
        test_mountain_wave()
        test_periodic_hill()
        test_vortex()
        test_wave_packet()
        test_wkb_mountain_wave()
        test_wkb_wave_packet()
    end

    @testset verbose = true "Code quality" begin
        test_formatting()

        test_ascii()
        test_julia_validity()
        test_literal_checks()
        test_long_named_function_definitions()
        test_type_checks()

        @testset "Aqua.jl" begin
            Aqua.test_all(PinCFlow)
        end
    end
end
