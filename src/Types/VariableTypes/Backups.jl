"""
```julia
Backups{A <: AbstractArray{<:AbstractFloat, 3}}
```

Container for backup copies needed in the semi-implicit time scheme.

```julia
Backups(domain::Domain)::Backups
```

Initialize backup arrays sized according to the dimensions of the MPI subdomain.

# Fields

  - `rhoold::A`: Density backup.

  - `rhopold::A`: Density-fluctuations backup.

  - `uold::A`: Zonal-wind backup.

  - `vold::A`: Meridional-wind backup.

  - `wold::A`: Transformed-vertical-wind backup.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Backups{A <: AbstractArray{<:AbstractFloat, 3}}
    rhoold::A
    rhopold::A
    uold::A
    vold::A
    wold::A
end

function Backups(domain::Domain)::Backups
    (; nxx, nyy, nzz) = domain

    return Backups([zeros(nxx, nyy, nzz) for i in 1:5]...)
end
