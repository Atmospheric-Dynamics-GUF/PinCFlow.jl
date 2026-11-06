"""
```julia
TestTools
```

Module that provides convenience functions for testing.

# Exports

  - [`TestTools.test_example`](@ref)
"""
module TestTools

using Test
using JuliaFormatter
using LinearAlgebra: norm
using MPI
using HDF5

include("compute_norms.jl")
include("replace_assignments.jl")
include("test_example.jl")

export test_example

end
