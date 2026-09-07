"""
```julia
Interpolation
```

Module for interpolating mean-flow quantities to ray-volume positions.

Provides functions that find the grid points closest to a given ray-volume position and perform trilinear interpolation of mean-flow quantities.

# See also

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Update`](@ref)
"""
module Interpolation

using ...Macros
using ...Types
using ...Update

"""
```julia
DUDX
```

Singleton for dispatch to interpolation of ``\\partial u_\\mathrm{b} / \\partial x``.
"""
struct DUDX end

"""
```julia
DUDY
```

Singleton for dispatch to interpolation of ``\\partial u_\\mathrm{b} / \\partial y``.
"""
struct DUDY end

"""
```julia
DUDZ
```

Singleton for dispatch to interpolation of ``\\partial u_\\mathrm{b} / \\partial z``.
"""
struct DUDZ end

"""
```julia
DVDX
```

Singleton for dispatch to interpolation of ``\\partial v_\\mathrm{b} / \\partial x``.
"""
struct DVDX end

"""
```julia
DVDY
```

Singleton for dispatch to interpolation of ``\\partial v_\\mathrm{b} / \\partial y``.
"""
struct DVDY end

"""
```julia
DVDZ
```

Singleton for dispatch to interpolation of ``\\partial v_\\mathrm{b} / \\partial z``.
"""
struct DVDZ end

"""
```julia
DX
```

Singleton for dispatch to interpolation of the gradient of a scalar field in ``\\hat{x}``-direction.
"""
struct DX end

"""
```julia
DY
```

Singleton for dispatch to interpolation of the gradient of a scalar field in ``\\hat{y}``-direction.
"""
struct DY end

"""
```julia
DZ
```

Singleton for dispatch to interpolation of the gradient of a scalar field in ``\\hat{z}``-direction.
"""
struct DZ end

include("compute_derivatives.jl")
include("get_next_half_level.jl")
include("get_next_level.jl")
include("interpolate_mean_flow.jl")
include("interpolate_scalar.jl")
include("interpolate.jl")

export DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DX, DY, DZ

export get_next_half_level,
    get_next_level, interpolate_mean_flow, interpolate_scalar

end
