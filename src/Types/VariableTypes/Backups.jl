"""
```julia
Backups{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
}
```

Container for backup copies needed in the semi-implicit time scheme.

```julia
Backups(domain::Domain)::Backups
```

Initialize backup arrays sized according to the dimensions of the MPI subdomain.

# Fields

  - `rhoold::A`: Density backup.

  - `rhopold::B`: Density-fluctuations backup.

  - `uold::C`: Zonal-wind backup.

  - `vold::D`: Meridional-wind backup.

  - `wold::E`: Transformed-vertical-wind backup.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Backups{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
}
    rhoold::A
    rhopold::B
    uold::C
    vold::D
    wold::E
end

function Backups(domain::Domain)::Backups
    (; nxx, nyy, nzz) = domain

    return Backups([zeros(nxx, nyy, nzz) for i in 1:5]...)
end
