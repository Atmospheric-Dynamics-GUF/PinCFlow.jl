"""
```julia
Rays{A <: AbstractArray{<:AbstractFloat, 4}}
```

Container for ray position, wavenumber, and propagation data.

# Fields

  - `x`, `y`, `z`: Ray positions in physical space
  - `k`, `l`, `m`: Ray wavenumbers in spectral space
  - `dxray`, `dyray`, `dzray`: Ray position increments
  - `dkray`, `dlray`, `dmray`: Ray wavenumber increments
  - `dens`: Ray density (wave action density)
"""
struct Rays{A <: AbstractArray{<:AbstractFloat, 4}}
    x::A
    y::A
    z::A
    k::A
    l::A
    m::A
    dxray::A
    dyray::A
    dzray::A
    dkray::A
    dlray::A
    dmray::A
    dens::A
end

"""
```julia
Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
```

Construct a `Rays` instance.

# Arguments

  - `nray_wrk`: Number of working rays
  - `nxx`: Number of grid points in x-direction
  - `nyy`: Number of grid points in y-direction
  - `nzz`: Number of grid points in z-direction

# Returns

  - `::Rays`: `Rays` instance with zero-initialized arrays.
"""
function Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
    return Rays([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:13]...)
end
