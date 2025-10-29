"""
```julia
WKBIncrements{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
    C <: AbstractArray{<:AbstractFloat, 4},
    D <: AbstractArray{<:AbstractFloat, 4},
    E <: AbstractArray{<:AbstractFloat, 4},
    F <: AbstractArray{<:AbstractFloat, 4},
    G <: AbstractArray{<:AbstractFloat, 4},
    H <: AbstractArray{<:AbstractFloat, 4},
    I <: AbstractArray{<:AbstractFloat, 4},
}
```

Ray-volume-propagation increments.

```julia
WKBIncrements(
    nray_wrk::Integer,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
)::WKBIncrements
```

Construct an `WKBIncrements` instance, with arrays sized according to the given dimensions.

# Fields

  - `dxray::A`: WKBIncrements for the position in ``x``.

  - `dyray::B`: WKBIncrements for the position in ``y``.

  - `dzray::C`: WKBIncrements for the position in ``z``.

  - `dkray::D`: WKBIncrements for the position in ``k``.

  - `dlray::E`: WKBIncrements for the position in ``l``.

  - `dmray::F`: WKBIncrements for the position in ``m``.

  - `ddxray::G`: WKBIncrements for the extent in ``x``.

  - `ddyray::H`: WKBIncrements for the extent in ``y``.

  - `ddzray::I`: WKBIncrements for the extent in ``z``.

# Arguments

  - `nray_wrk`: Size of the spectral dimension of ray-volume arrays.

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct WKBIncrements{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
    C <: AbstractArray{<:AbstractFloat, 4},
    D <: AbstractArray{<:AbstractFloat, 4},
    E <: AbstractArray{<:AbstractFloat, 4},
    F <: AbstractArray{<:AbstractFloat, 4},
    G <: AbstractArray{<:AbstractFloat, 4},
    H <: AbstractArray{<:AbstractFloat, 4},
    I <: AbstractArray{<:AbstractFloat, 4},
}
    dxray::A
    dyray::B
    dzray::C
    dkray::D
    dlray::E
    dmray::F
    ddxray::G
    ddyray::H
    ddzray::I
end

function WKBIncrements(
    nray_wrk::Integer,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
)::WKBIncrements
    return WKBIncrements([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:9]...)
end
