"""
```julia
RayUpdate
```

Module for the integration of the ray equations.

In addition to ray-volume initialization and propagation, functions for tracking ray volumes on the model grid and controlling their count, as well as a saturation scheme for capturing wave breaking, are provided.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations`](@ref)

  - [`PinCFlow.MSGWaM.RaySources`](@ref)
"""
module RayUpdate

using MPI
using Statistics
using ..BoundaryRays
using ..Interpolation
using ..RayOperations
using ..RaySources
using ...Types
using ...PinCFlow

"""
```julia
X
```

Singleton for dispatch to operations in ``x``-direction.
"""
struct X end

"""
```julia
Y
```

Singleton for dispatch to operations in ``y``-direction.
"""
struct Y end

"""
```julia
Z
```

Singleton for dispatch to operations in ``z``-direction.
"""
struct Z end

"""
```julia
XZ
```

Singleton for dispatch to operations in ``x``- and ``z``-direction.
"""
struct XZ end

"""
```julia
YZ
```

Singleton for dispatch to operations in ``y``- and ``z``-direction.
"""
struct YZ end

"""
```julia
XYZ
```

Singleton for dispatch to operations in all directions.
"""
struct XYZ end

include("apply_saturation_scheme!.jl")
include("initialize_rays!.jl")
include("merge_rays!.jl")
include("propagate_rays!.jl")
include("shift_rays!.jl")
include("split_rays!.jl")

export X, Y, Z, XZ, YZ, XYZ

export apply_saturation_scheme!,
    initialize_rays!, merge_rays!, propagate_rays!, shift_rays!, split_rays!

end
