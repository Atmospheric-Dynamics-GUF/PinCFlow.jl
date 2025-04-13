module PoissonTypes

using ..FoundationalTypes

include("Tensor.jl")
include("Operator.jl")
include("Preconditioner.jl")
include("BicGStab.jl")
include("Correction.jl")
include("Poisson.jl")

export Tensor, Operator, Preconditioner, BicGStab, Correction, Poisson

end