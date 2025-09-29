"""
```julia
set_plot_style()
```

Configure PythonPlot.jl to use a preset plot style.
"""
function set_plot_style end

"""
```julia
symmetric_contours(
    minimum::AbstractFloat,
    maximum::AbstractFloat;
    number::Integer = 10,
    colormap_name::String = "seismic",
)::Tuple{<:LinRange{<:AbstractFloat, <:Integer}, <:Any}
```

Compute symmetric contours levels and return them and a correspondingly cropped colormap.

# Arguments

  - `minimum`: Smallest value to be plotted.

  - `maximum`: Largest value to be plotted.

# Keywords

  - `number`: Number of contour levels.

  - `colormap_name`: Name under which the chosen colormap is registered.
"""
function symmetric_contours end