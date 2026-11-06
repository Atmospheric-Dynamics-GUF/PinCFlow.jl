using CairoMakie
using PinCFlow
using Test

using .Meta: parse, parseall
using MacroTools: postwalk
using REPL: moduleusings

const PinCFlowMakieExt = Base.get_extension(PinCFlow, :PinCFlowMakieExt)

include("../test/test_docstrings.jl")
include("../test/test_files.jl")

test_docstrings(; fix = true, folder = "src")
test_docstrings(; fix = true, folder = "ext")
