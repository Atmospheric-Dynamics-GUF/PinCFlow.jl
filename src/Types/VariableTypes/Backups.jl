"""
```julia
Backups{A <: AbstractArray{<:AbstractFloat, 3}}
```

Container for backup copies needed in the semi-implicit time scheme.

```julia
Backups(predictands::Predictands)::Backups
```

Initialize backup arrays sized according to the dimensions of the arrays in `predictands`.

# Fields

  - `rhoold::A`: Density backup.

  - `rhopold::A`: Density-fluctuations backup.

  - `uold::A`: Zonal-wind backup.

  - `vold::A`: Meridional-wind backup.

  - `wold::A`: Transformed-vertical-wind backup.

  - `pipold::A`: Exner-pressure-fluctuations backup.

  - `p::A`: Mass-weighted-potential-temperature backup.

# Arguments

  - `predictands`: Predictands for which the backups are needed.
"""
struct Backups{A <: AbstractArray{<:AbstractFloat, 3}}
    rhoold::A
    rhopold::A
    uold::A
    vold::A
    wold::A
    pipold::A
    pold::A
end

function Backups(predictands::Predictands)::Backups
    return Backups(
        (
            zeros(size(getfield(predictands, field))) for
            field in fieldnames(Predictands)
        )...,
    )
end
