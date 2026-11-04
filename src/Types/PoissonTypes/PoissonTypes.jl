"""
```julia
PoissonTypes
```

Module for composite types used by the Poisson solver.

# Dependencies

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)

# Exports

  - [`PinCFlow.Types.PoissonTypes.BiCGSTAB`](@ref)

  - [`PinCFlow.Types.PoissonTypes.Correction`](@ref)

  - [`PinCFlow.Types.PoissonTypes.Operator`](@ref)

  - [`PinCFlow.Types.PoissonTypes.Poisson`](@ref)

  - [`PinCFlow.Types.PoissonTypes.Tensor`](@ref)
"""
module PoissonTypes

using ..NamelistTypes
using ..FoundationalTypes
using ...Macros

include("Tensor.jl")
include("Operator.jl")
include("BiCGSTAB.jl")
include("Correction.jl")
include("Poisson.jl")

export Tensor, Operator, BiCGSTAB, Correction, Poisson

end
