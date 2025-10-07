"""
```julia
TestTools
```

Extension that provides convenience functions for testing.
"""
module TestTools

using Test
using LinearAlgebra: norm
using TrixiTest: trixi_include, get_kwarg
using MPI
using HDF5

include("@test_example.jl")
include("compute_norms.jl")
include("test_example.jl")

export @test_example, test_example

end
