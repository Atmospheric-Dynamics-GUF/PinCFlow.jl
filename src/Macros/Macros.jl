"""
```julia
Macros
```

Module that contains PinCFlow.jl's convenience macros.
"""
module Macros

import Polyester
import FastBroadcast

include("dispatch.jl")
include("find_argument.jl")
include("ivy.jl")
include("replace_argument.jl")

include("@dispatch.jl")
include("@ivy.jl")
include("@last.jl")
include("@share.jl")

export @dispatch, @ivy, @share

end
