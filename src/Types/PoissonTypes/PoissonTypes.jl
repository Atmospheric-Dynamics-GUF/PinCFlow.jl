"""
```julia
PoissonTypes
```

Module for composite types used by the Poisson solver.
"""
module PoissonTypes

using ..NamelistTypes
using ..FoundationalTypes

include("Tensor.jl")
include("Operator.jl")
include("Preconditioner.jl")
include("BicGStab.jl")
include("Correction.jl")
include("Poisson.jl")

export Tensor, Operator, Preconditioner, BicGStab, Correction, Poisson

end
