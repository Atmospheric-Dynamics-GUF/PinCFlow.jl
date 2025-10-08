"""
```julia
GridNamelist{A <: AbstractFloat, B <: Integer}
```

Namelist for parameters describing the grid.

```julia
GridNamelist(;
    mountain_height::AbstractFloat = 1.0E+2,
    mountain_half_width::AbstractFloat = 1.0E+3,
    mountain_case::Integer = 1,
    height_factor::AbstractFloat = 1.0E+0,
    width_factor::AbstractFloat = 1.0E+0,
    spectral_modes::Integer = 1,
    stretch_exponent::AbstractFloat = 1.0E+0,
)::GridNamelist
```

Construct a `GridNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `mountain_height::A`: Mountain height (summit).

  - `mountain_half_width::A`: Mountain half-width.

  - `mountain_case::B`: Mountain shape. See [`PinCFlow.Types.FoundationalTypes.compute_topography`](@ref) for available options.

  - `height_factor::A`: Factor between the amplitudes of the large and small orographic features.

  - `width_factor::A`: Factor between the wavelengths of the large and small orographic features.

  - `spectral_modes::B`: Number of spectral modes used for `mountain_case == 13`.

  - `stretch_exponent::A`: Vertical-grid-stretching parameter.
"""
struct GridNamelist{A <: AbstractFloat, B <: Integer}
    mountain_height::A
    mountain_half_width::A
    mountain_case::B
    height_factor::A
    width_factor::A
    spectral_modes::B
    stretch_exponent::A
end

function GridNamelist(;
    mountain_height::AbstractFloat = 1.0E+2,
    mountain_half_width::AbstractFloat = 1.0E+3,
    mountain_case::Integer = 1,
    height_factor::AbstractFloat = 1.0E+0,
    width_factor::AbstractFloat = 1.0E+0,
    spectral_modes::Integer = 1,
    stretch_exponent::AbstractFloat = 1.0E+0,
)::GridNamelist
    return GridNamelist(
        mountain_height,
        mountain_half_width,
        mountain_case,
        height_factor,
        width_factor,
        spectral_modes,
        stretch_exponent,
    )
end
