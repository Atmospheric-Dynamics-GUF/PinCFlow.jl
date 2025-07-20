"""
```julia
GridNamelist{A <: AbstractFloat, B <: Integer}
```

Namelist for parameters describing the grid.

# Fields

  - `mountainheight_dim::A`: Mountain height (summit).
  - `mountainwidth_dim::A`: Mountain half-width.
  - `mountain_case::B`: Mountain shape. See [`PinCFlow.Types.FoundationalTypes.compute_topography`](@ref) for available options.
  - `height_factor::A`: Factor between the amplitudes of the large and small orographic features.
  - `width_factor::A`: Factor between the wavelengths of the large and small orographic features.
  - `spectral_modes::B`: Number of spectral modes used for `mountain_case == 13`.
  - `stretch_exponent::A`: Vertical-grid-stretching parameter.
"""
struct GridNamelist{A <: AbstractFloat, B <: Integer}
    mountainheight_dim::A
    mountainwidth_dim::A
    mountain_case::B
    height_factor::A
    width_factor::A
    spectral_modes::B
    stretch_exponent::A
end

"""
```julia
GridNamelist(;
    mountainheight_dim::AbstractFloat = 1.0E+2,
    mountainwidth_dim::AbstractFloat = 1.0E+3,
    mountain_case::Integer = 1,
    height_factor::AbstractFloat = 1.0E+0,
    width_factor::AbstractFloat = 1.0E+0,
    spectral_modes::Integer = 1,
    stretch_exponent::AbstractFloat = 1.0E+0,
)
```

Construct a `GridNamelist` instance with the given keyword arguments as properties.

# Arguments

  - `mountainheight_dim`: Mountain height (summit).
  - `mountainwidth_dim`: Mountain half-width.
  - `mountain_case`: Mountain shape. See [`PinCFlow.Types.FoundationalTypes.compute_topography`](@ref) for available options.
  - `height_factor`: Factor between the amplitudes of the large and small orographic features.
  - `width_factor`: Factor between the wavelengths of the large and small orographic features.
  - `spectral_modes`: Number of spectral modes used for `mountain_case == 13`.
  - `stretch_exponent`: Vertical-grid-stretching parameter.

# Returns

  - `::GridNamelist`: `GridNamelist` instance.
"""
function GridNamelist(;
    mountainheight_dim::AbstractFloat = 1.0E+2,
    mountainwidth_dim::AbstractFloat = 1.0E+3,
    mountain_case::Integer = 1,
    height_factor::AbstractFloat = 1.0E+0,
    width_factor::AbstractFloat = 1.0E+0,
    spectral_modes::Integer = 1,
    stretch_exponent::AbstractFloat = 1.0E+0,
)
    return GridNamelist(
        mountainheight_dim,
        mountainwidth_dim,
        mountain_case,
        height_factor,
        width_factor,
        spectral_modes,
        stretch_exponent,
    )
end
