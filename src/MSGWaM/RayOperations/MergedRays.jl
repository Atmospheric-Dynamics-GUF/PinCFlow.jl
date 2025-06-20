"""
    MergedRays{A, B}

# Type Parameters

  - `A <: AbstractVector{<:AbstractFloat}`: Vector type for position and extent arrays
  - `B <: Ref{<:AbstractFloat}`: Reference type for scalar values

# Fields

  - `xr::A`: X-direction position bounds
  - `yr::A`: Y-direction position bounds
  - `zr::A`: Z-direction position bounds
  - `kr::A`: K-direction wavenumber bounds
  - `lr::A`: L-direction wavenumber bounds
  - `mr::A`: M-direction wavenumber bounds
  - `nr::B`: Wave action density reference
"""
struct MergedRays{
    A <: AbstractVector{<:AbstractFloat},
    B <: Ref{<:AbstractFloat},
}
    xr::A
    yr::A
    zr::A
    kr::A
    lr::A
    mr::A
    nr::B
end
