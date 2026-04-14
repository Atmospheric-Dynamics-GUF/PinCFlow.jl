module Smoothing

using ...Types
using ...Boundaries
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

include("smooth_gw_tendencies!.jl")
include("smooth_gw_amplitudes!.jl")
include("apply_shapiro_filter!.jl")

export X, Y, Z, XZ, YZ, XYZ

export smooth_gw_amplitudes!, smooth_gw_tendencies!

end