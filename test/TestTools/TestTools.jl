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
using MPI
using HDF5

include("compute_norms.jl")
include("occurrences.jl")
include("replace_assignments.jl")
include("test_example.jl")

export occurrences, test_example

end
