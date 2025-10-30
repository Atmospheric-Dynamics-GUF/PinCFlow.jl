"""
```julia
PoissonTypes
```

Module for composite types used by the Poisson solver.

# See also

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)
"""
module PoissonTypes

using ..NamelistTypes
using ..FoundationalTypes
using ...PinCFlow

include("Tensor.jl")
include("Operator.jl")
include("BicGStab.jl")
include("Correction.jl")
include("Poisson.jl")

export Tensor, Operator, BicGStab, Correction, Poisson

end
