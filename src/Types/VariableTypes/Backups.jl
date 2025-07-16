"""
```julia
Backups{A <: AbstractArray{<:AbstractFloat, 3}}
```

Storage container for backup copies of predictand variables during time stepping.

# Fields

  - `rhoold::A`: Previous density values (nxx x nyy x nzz)
  - `rhopold::A`: Previous density fluctuation values (nxx x nyy x nzz)
  - `uold::A`: Previous zonal velocity values (nxx x nyy x nzz)
  - `vold::A`: Previous meridional velocity values (nxx x nyy x nzz)
  - `wold::A`: Previous vertical velocity values (nxx x nyy x nzz)

# Usage

Backup fields store previous time step values for multi-stage time integration schemes.
Used by [`PinCFlow.Integration.save_backups!`](@ref) to preserve state before updates.
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

Initialize backup storage arrays sized according to domain decomposition.

# Arguments

  - `domain`: Domain specification containing grid dimensions

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
