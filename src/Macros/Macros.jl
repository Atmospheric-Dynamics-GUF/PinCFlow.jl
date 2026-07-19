"""
```julia
Macros
```

Module that contains PinCFlow.jl's convenience macros.
"""
module Macros

include("find_argument.jl")
include("ivy.jl")
include("replace_argument.jl")

include("@dispatch.jl")
include("@ivy.jl")

export @dispatch, @ivy

end
