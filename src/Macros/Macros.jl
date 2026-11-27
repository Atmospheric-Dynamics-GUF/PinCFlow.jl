"""
```julia
Macros
```

Module that contains PinCFlow.jl's convenience macros.
"""
module Macros

import FastBroadcast
using MacroTools: isdef, longdef

include("find_argument.jl")
include("ivy.jl")
include("replace_argument.jl")

include("@dispatch.jl")
include("@fuse.jl")
include("@ivy.jl")

export @dispatch, @fuse, @ivy

end
