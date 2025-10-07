"""
```julia
plot_contours(
    file::AbstractString,
    data::HDF5.File,
    variable::AbstractString,
    label::AbstractString,
    indices::Vararg{<:NTuple{4, <:Integer}};
    number::Integer = 10,
    colormap_name::Symbol = :seismic,
)
```

Create contour plots of the dataset `variable` in `data`, display it and save it to `file`.

# Arguments

  - `file`: File to save the plots to.

  - `data`: PinCFlow output data.

  - `variable`: Output variable to plot.

  - `label`: Colorbar label for the plots.

  - `indices`: Tuples of four indices. The first three indices of each tuple define the planes in which the contours are to be plotted, whereas the fourth is the temporal index.

# Keywords

  - `number`: Number of contour levels.

  - `colormap_name`: Colormap of choice.
"""
function plot_contours end
