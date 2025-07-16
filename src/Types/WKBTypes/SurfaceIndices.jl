"""
```julia
SurfaceIndices{A <: AbstractArray{<:Integer, 3}, B <: AbstractVector{<:Integer}}
```

Indices for ray launching at surface boundaries.

# Fields

  - `ir_sfc`: Surface ray indices (3D array)
  - `ix2_sfc`, `jy2_sfc`, `kz2_sfc`: Grid position indices for surface rays
  - `ik_sfc`, `jl_sfc`, `km_sfc`: Wavenumber indices for surface rays
  - `iwm_sfc`: Wave mode indices for surface rays
"""
struct SurfaceIndices{
    A <: AbstractArray{<:Integer, 3},
    B <: AbstractVector{<:Integer},
}
    ir_sfc::A
    ix2_sfc::B
    jy2_sfc::B
    kz2_sfc::B
    ik_sfc::B
    jl_sfc::B
    km_sfc::B
    iwm_sfc::B
end

"""
```julia
SurfaceIndices(n_sfc::Integer, nxx::Integer, nyy::Integer)
```

Construct a `SurfaceIndices` instance.

# Arguments

  - `n_sfc`: Number of surface rays
  - `nxx`: Number of grid points in x-direction
  - `nyy`: Number of grid points in y-direction

# Returns

  - `::SurfaceIndices`: `SurfaceIndices` instance with zero-initialized arrays.
"""
function SurfaceIndices(n_sfc::Integer, nxx::Integer, nyy::Integer)
    return SurfaceIndices(
        zeros(Int, n_sfc, nxx, nyy),
        [zeros(Int, n_sfc) for i in 1:7]...,
    )
end
