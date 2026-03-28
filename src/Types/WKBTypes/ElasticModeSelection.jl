"""
```julia
struct ElasticModeSelection{
    A <: AbstractVector{<:Integer},
    B <: AbstractMatrix{<:Integer},
    C <: AbstractMatrix{<:AbstractFloat},
}
```

Composite type for elastic-mode-selection data.

```julia
ElasticModeSelection(
    wave_modes::Integer,
    nxx::Integer,
    nyy::Integer,
)::ElasticModeSelection
```

Construct an `ElasticModeSelection` instance with arrays sized according to the given dimensions.

# Fields

  - `sorted_wave_mode_indices::A`: Array for indices that sort the wave modes.

  - `launch_mode_count::B`: Array that stores the numbers of selected modes.

  - `launch_power_fraction::C`: Array that stores the power fractions retained by the selection algorithm.
"""
struct ElasticModeSelection{
    A <: AbstractVector{<:Integer},
    B <: AbstractMatrix{<:Integer},
    C <: AbstractMatrix{<:AbstractFloat},
}
    sorted_wave_mode_indices::A
    launch_mode_count::B
    launch_power_fraction::C
end

function ElasticModeSelection(
    wave_modes::Integer,
    nxx::Integer,
    nyy::Integer,
)::ElasticModeSelection
    return ElasticModeSelection(
        zeros(Int, wave_modes),
        zeros(Int, nxx, nyy),
        zeros(nxx, nyy),
    )
end
