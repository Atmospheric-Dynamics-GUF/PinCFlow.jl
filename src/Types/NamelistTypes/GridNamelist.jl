"""
```julia
GridNamelist{A <: AbstractFloat, B <: Integer, C <: Function, D <: Function}
```

Namelist for parameters describing the grid.

```julia
GridNamelist(;
    mountain_height::AbstractFloat = 1.0E+2,
    mountain_half_width::AbstractFloat = 1.0E+3,
    height_factor::AbstractFloat = 1.0E+0,
    width_factor::AbstractFloat = 1.0E+0,
    spectral_modes::Integer = 1,
    stretch_exponent::AbstractFloat = 1.0E+0,
    resolved_topography::Function = mountain_range,
    unresolved_topography::Function = mountain_range,
)::GridNamelist
```

Construct a `GridNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `mountain_height::A`: Mountain height (summit) used by `mountain_range`.

  - `mountain_half_width::A`: Mountain half-width used by `mountain_range`.

  - `height_factor::A`: Factor between the amplitudes of the large and small orographic features in `mountain_range`.

  - `width_factor::A`: Factor between the wavelengths of the large and small orographic features in `mountain_range`.

  - `spectral_modes::B`: Number of spectral modes used by `mountain_range`.

  - `stretch_exponent::A`: Vertical-grid-stretching parameter.

  - `resolved_topography::C`: Function that returns the resolved topography at a specified horizontal position (see [`PinCFlow.Types.NamelistTypes.mountain_range`](@ref) for the default option).

  - `resolved_topography::D`: Function that returns a specified spectral mode of the unresolved topography at a specified horizontal position (see [`PinCFlow.Types.NamelistTypes.mountain_range`](@ref) for the default option).
"""
struct GridNamelist{
    A <: AbstractFloat,
    B <: Integer,
    C <: Function,
    D <: Function,
}
    mountain_height::A
    mountain_half_width::A
    height_factor::A
    width_factor::A
    spectral_modes::B
    stretch_exponent::A
    resolved_topography::C
    unresolved_topography::D
end

function GridNamelist(;
    mountain_height::AbstractFloat = 1.0E+2,
    mountain_half_width::AbstractFloat = 1.0E+3,
    height_factor::AbstractFloat = 1.0E+0,
    width_factor::AbstractFloat = 1.0E+0,
    spectral_modes::Integer = 1,
    stretch_exponent::AbstractFloat = 1.0E+0,
    resolved_topography::Function = mountain_range,
    unresolved_topography::Function = mountain_range,
)::GridNamelist
    return GridNamelist(
        mountain_height,
        mountain_half_width,
        height_factor,
        width_factor,
        spectral_modes,
        stretch_exponent,
        resolved_topography,
        unresolved_topography,
    )
end
