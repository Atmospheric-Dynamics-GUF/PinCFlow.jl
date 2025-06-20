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
GridNamelist(; <keyword arguments>)

Configuration parameters for topography and vertical grid stretching.

# Arguments

  - `mountainheight_dim::AbstractFloat = 1.0E+2`: Mountain height [m]
  - `mountainwidth_dim::AbstractFloat = 1.0E+3`: Mountain width [m]
  - `mountain_case::Integer = 1`: Mountain profile selector. See [compute_topography] for a available options.
  - `height_factor::AbstractFloat = 1.0E+0`: Topography height scaling
  - `width_factor::AbstractFloat = 1.0E+0`: Topography width scaling
  - `spectral_modes::Integer = 1`: Number of spectral modes for terrain representation
  - `stretch_exponent::AbstractFloat = 1.0E+0`: Vertical grid stretching parameter
"""
function GridNamelist(;
    mountainheight_dim = 1.0E+2,
    mountainwidth_dim = 1.0E+3,
    mountain_case = 1,
    height_factor = 1.0E+0,
    width_factor = 1.0E+0,
    spectral_modes = 1,
    stretch_exponent = 1.0E+0,
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
