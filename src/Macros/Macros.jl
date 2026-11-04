"""
```julia
Macros
```

Module that contains PinCFlow.jl's convenience macros.

# Exports

  - [`PinCFlow.Macros.@dispatch`](@ref)

  - [`PinCFlow.Macros.@ivy`](@ref)
"""
module Macros

include("find_argument.jl")
include("ivy.jl")
include("replace_argument.jl")

include("@dispatch.jl")
include("@ivy.jl")

export @dispatch, @ivy

end
