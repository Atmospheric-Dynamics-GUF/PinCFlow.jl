"""
```julia
symmetric_contours(
    minimum::AbstractFloat,
    maximum::AbstractFloat;
    number::Integer = 10,
    colormap_name::Symbol = :seismic,
)::Tuple{<:LinRange{<:AbstractFloat, <:Integer}, <:Any}
```

Compute symmetric contours levels and return them and a correspondingly indexed colormap.

# Arguments

  - `minimum`: Smallest value to be plotted.

  - `maximum`: Largest value to be plotted.

# Keywords

  - `number`: Number of contour levels.

  - `colormap_name`: Name under which the chosen colormap is registered.
"""
function symmetric_contours end
