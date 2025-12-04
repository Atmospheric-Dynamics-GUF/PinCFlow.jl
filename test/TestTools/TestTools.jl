"""
```julia
TestTools
```

Module that provides convenience functions for testing.
"""
module TestTools

using Test
using JuliaFormatter
using LinearAlgebra: norm
using PinCFlow: replace_assignments
using MPI
using HDF5

include("compute_norms.jl")
include("test_example.jl")

export test_example

end
