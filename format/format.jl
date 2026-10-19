using JuliaFormatter
using Test

include("../test/test_formatting.jl")

test_formatting(; directory = ".", fix = true)
