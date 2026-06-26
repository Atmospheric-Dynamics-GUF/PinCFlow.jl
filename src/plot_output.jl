"""
```julia
plot_output(
    plot_file::AbstractString,
    data_file::AbstractString,
    fields::Vararg{
        Union{
            Tuple{Symbol, <:Integer},
            Tuple{Symbol, <:Real, <:Real, <:Real, <:Integer},
        },
    };
    colormap_name::Symbol = :seismic,
    number::Integer = 10,
    significant_digits::Integer = 3,
    space_unit::Symbol = :km,
    time_unit::Symbol = :h,
)
```

Create contour plots of the dataset `variable` in `data`, display it and save it to `file`.

# Arguments

  - `plot_file`: File to save the plots to.

  - `data_file`: HDF5 file with PinCFlow.jl output data.

  - `fields`: Either tuples of a variable name and a temporal index (for 2D data), or tuples of a variable name, three fractions which define the relative positions of the ``\\hat{y}``-``\\hat{z}``, ``\\hat{x}``-``\\hat{z}``, and ``\\hat{x}``-``\\hat{y}`` planes, and a temporal index (for 3D data).

# Keywords

  - `colormap_name`: Colormap of choice.

  - `number`: Number of contour levels.

  - `significant_digits`: Significant digits to which the output data is rounded before visualization.

  - `space_unit`: Unit used for the coordinates. Must be `:km` or `:m`.

  - `time_unit`: Unit used for the time. Must be `:d`, `:h`, `:min` or `:s`.
"""
function plot_output end
