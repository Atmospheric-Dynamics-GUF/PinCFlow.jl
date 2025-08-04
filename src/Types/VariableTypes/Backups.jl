"""
```julia
Backups{A <: AbstractArray{<:AbstractFloat, 3}}
```

Container for backup copies needed in the semi-implicit time scheme.

# Fields

  - `rhoold::A`: Density backup.
  - `rhopold::A`: Density-fluctuations backup.
  - `uold::A`: Zonal-wind backup.
  - `vold::A`: Meridional-wind backup.
  - `wold::A`: Transformed-vertical-wind backup.
"""
struct Backups{A <: AbstractArray{<:AbstractFloat, 3}}
    rhoold::A
    rhopold::A
    uold::A
    vold::A
    wold::A
end

"""
```julia
Backups(domain::Domain)
```

Initialize backup arrays sized according to the dimensions of the MPI subdomain.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Returns

  - `::Backups`: `Backups` instance with zero-initialized arrays.
"""
function Backups(domain::Domain)
    (; nxx, nyy, nzz) = domain

    # Initialize the backups.
    (rhoold, rhopold, uold, vold, wold) = (zeros(nxx, nyy, nzz) for i in 1:5)

    # Return a Backups instance.
    return Backups(rhoold, rhopold, uold, vold, wold)
end
