"""
```julia
ElasticModeSelection{
    A <: AbstractVector{<:Integer},
    B <: AbstractMatrix{<:Integer},
    C <: AbstractMatrix{<:AbstractFloat},
}
```

Composite type for elastic-mode-selection data.

```julia
ElasticModeSelection(
    integer_type::DataType,
    float_type::DataType,
    wave_modes::Integer,
    nxx::Integer,
    nyy::Integer,
)::ElasticModeSelection
```

Construct an `ElasticModeSelection` instance with arrays sized according to the given dimensions.

# Fields

  - `integer_type`: Data type of `sorted_wave_mode_indices` and `launch_mode_count`.

  - `float_type`: Data type of `launch_power_fraction`.

  - `sorted_wave_mode_indices::A`: Array for indices that sort the wave modes.

  - `launch_mode_count::B`: Array that stores the numbers of selected modes.

  - `launch_power_fraction::C`: Array that stores the power fractions retained by the selection algorithm.

!!! danger "Experimental"
    The elastic mode selection is an experimental feature adapted from [Banerjee (2026)](https://doi.org/10.5281/zenodo.20582010).
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
    integer_type::DataType,
    float_type::DataType,
    wave_modes::Integer,
    nxx::Integer,
    nyy::Integer,
)::ElasticModeSelection
    return ElasticModeSelection(
        zeros(integer_type, wave_modes),
        zeros(integer_type, nxx, nyy),
        zeros(float_type, nxx, nyy),
    )
end
