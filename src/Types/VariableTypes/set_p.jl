"""
```julia
set_p(
    model::AbstractModel,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
```

Return a zero-size array in non-compressible modes.

In these cases, the mass-weighted potential temperature is a background field: constant in Boussinesq mode, vertically varying in pseudo-incompressible mode.

```julia
set_p(
    model::Compressible,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
```

Return a copy of `pbar` in compressible mode.

In compressible mode, the mass-weighted potential temperature is a prognostic variable.

# Arguments:

  - `mode`: Dynamic equations.

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `pbar`: Mass-weighted potential temperature.

"""
function set_p end

function set_p(
    model::AbstractModel,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return zeros(0, 0, 0)
end

function set_p(
    model::Compressible,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return copy(pbar)
end
