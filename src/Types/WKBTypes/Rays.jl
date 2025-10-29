"""
```julia
Rays{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
    C <: AbstractArray{<:AbstractFloat, 4},
    D <: AbstractArray{<:AbstractFloat, 4},
    E <: AbstractArray{<:AbstractFloat, 4},
    F <: AbstractArray{<:AbstractFloat, 4},
    G <: AbstractArray{<:AbstractFloat, 4},
    H <: AbstractArray{<:AbstractFloat, 4},
    I <: AbstractArray{<:AbstractFloat, 4},
    J <: AbstractArray{<:AbstractFloat, 4},
    K <: AbstractArray{<:AbstractFloat, 4},
    L <: AbstractArray{<:AbstractFloat, 4},
    M <: AbstractArray{<:AbstractFloat, 4},
}
```

Container for prognostic ray-volume properties.

```julia
Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)::Rays
```

Construct a `Rays` instance, with arrays sized according to the given dimensions.

# Fields

  - `x::A`: Position in ``x``.

  - `y::B`: Position in ``y``.

  - `z::C`: Position in ``z``.

  - `k::D`: Position in ``k``.

  - `l::E`: Position in ``l``.

  - `m::F`: Position in ``m``.

  - `dxray::G`: Extent in ``x``.

  - `dyray::H`: Extent in ``y``.

  - `dzray::I`: Extent in ``z``.

  - `dkray::J`: Extent in ``k``.

  - `dlray::K`: Extent in ``l``.

  - `dmray::L`: Extent in ``m``.

  - `dens::M`: Phase-space wave-action density.

# Arguments

  - `nray_wrk`: Size of the spectral dimension of ray-volume arrays.

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct Rays{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
    C <: AbstractArray{<:AbstractFloat, 4},
    D <: AbstractArray{<:AbstractFloat, 4},
    E <: AbstractArray{<:AbstractFloat, 4},
    F <: AbstractArray{<:AbstractFloat, 4},
    G <: AbstractArray{<:AbstractFloat, 4},
    H <: AbstractArray{<:AbstractFloat, 4},
    I <: AbstractArray{<:AbstractFloat, 4},
    J <: AbstractArray{<:AbstractFloat, 4},
    K <: AbstractArray{<:AbstractFloat, 4},
    L <: AbstractArray{<:AbstractFloat, 4},
    M <: AbstractArray{<:AbstractFloat, 4},
}
    x::A
    y::B
    z::C
    k::D
    l::E
    m::F
    dxray::G
    dyray::H
    dzray::I
    dkray::J
    dlray::K
    dmray::L
    dens::M
end

function Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)::Rays
    return Rays([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:13]...)
end
