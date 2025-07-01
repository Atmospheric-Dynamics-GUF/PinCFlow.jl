"""
```julia
Interpolation
```

Module for interpolating mean-flow quantities to ray-volume positions.

Provides functions that find the grid points closest to a given ray-volume position and perform trilinear interpolation of mean-flow quantities.
"""
module Interpolation

using ...Types
using ...Update

"""
```julia
N2 <: AbstractVariable
```

Singleton for dispatch to interpolation of ``𝑁²``.
"""
struct N2 <: AbstractVariable end

"""
```julia
DN2DZ <: AbstractVariable
```

Singleton for dispatch to interpolation of ``∂𝑁²/∂z``.
"""
struct DN2DZ <: AbstractVariable end

"""
```julia
DUDX <: AbstractVariable
```

Singleton for dispatch to interpolation of ``∂𝑢/∂𝑥``.
"""
struct DUDX <: AbstractVariable end

"""
```julia
DUDY <: AbstractVariable
```

Singleton for dispatch to interpolation of ``∂𝑢/∂𝑦``.
"""
struct DUDY <: AbstractVariable end

"""
```julia
DUDZ <: AbstractVariable
```

Singleton for dispatch to interpolation of ``∂𝑢/∂𝑧``.
"""
struct DUDZ <: AbstractVariable end

"""
```julia
DVDX <: AbstractVariable
```

Singleton for dispatch to interpolation of ``∂𝑣/∂𝑥``.
"""
struct DVDX <: AbstractVariable end

"""
```julia
DVDY <: AbstractVariable
```

Singleton for dispatch to interpolation of ``∂𝑣/∂𝑦``.
"""
struct DVDY <: AbstractVariable end

"""
```julia
DVDZ <: AbstractVariable
```

Singleton for dispatch to interpolation of ``∂𝑣/∂𝑧``.
"""
struct DVDZ <: AbstractVariable end

include("compute_derivatives.jl")
include("get_next_half_level.jl")
include("get_next_level.jl")
include("interpolate_mean_flow.jl")
include("interpolate_sponge.jl")
include("interpolate_stratification.jl")
include("interpolate.jl")

export N2, DN2DZ, DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ

export get_next_half_level,
    get_next_level,
    interpolate_mean_flow,
    interpolate_sponge,
    interpolate_stratification

end
