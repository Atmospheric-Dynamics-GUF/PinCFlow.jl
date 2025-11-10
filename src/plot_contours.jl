"""
```julia
plot_contours(
    file::AbstractString,
    data::HDF5.File,
    variable::AbstractString,
    indices::Vararg{NTuple{4, <:Integer}};
    number::Integer = 10,
    colormap_name::Symbol = :seismic,
    label::AbstractString = "",
    space_unit::AbstractString = "km",
    time_unit::AbstractString = "h",
)
```

Create contour plots of the dataset `variable` in `data`, display it and save it to `file`.

# Arguments

  - `file`: File to save the plots to.

  - `data`: PinCFlow.jl output data.

  - `variable`: Output variable to plot.

  - `indices`: Tuples of four indices. The first three indices of each tuple define the planes in which the contours are to be plotted, whereas the fourth is the temporal index.

# Keywords

  - `number`: Number of contour levels.

  - `colormap_name`: Colormap of choice.

  - `label`: Colorbar label for the plots.

  - `space_unit`: Unit used for the coordinates. Must be `"km"` or `"m"`.

  - `time_unit`: Unit used for the time. Must be `"d"`, `"h"`, `"min"` or `"s"`.
"""
function plot_contours end
