"""
```julia
plot_output(
    file::AbstractString,
    data::HDF5.File,
    fields::Vararg{
        Tuple{<:AbstractString, <:Integer, <:Integer, <:Integer, <:Integer},
    };
    colormap_name::Symbol = :seismic,
    number::Integer = 10,
    space_unit::Symbol = :km,
    time_unit::Symbol = :h,
)
```

Create contour plots of the dataset `variable` in `data`, display it and save it to `file`.

# Arguments

  - `file`: File to save the plots to.

  - `data`: PinCFlow.jl output data.

  - `fields`: Tuples of a variable name and four indices. The first three indices of each tuple define the planes in which the contours are to be plotted, whereas the fourth is the temporal index.

# Keywords

  - `number`: Number of contour levels.

  - `colormap_name`: Colormap of choice.

  - `label`: Colorbar label for the plots.

  - `space_unit`: Unit used for the coordinates. Must be `:km` or `:m`.

  - `time_unit`: Unit used for the time. Must be `:d`, `:h`, `:min` or `:s`.
"""
function plot_output end
