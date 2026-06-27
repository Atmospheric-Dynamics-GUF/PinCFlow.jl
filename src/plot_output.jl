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
    color_tick_format::AbstractString = "{:.1E}",
    number::Integer = 10,
    space_unit::Symbol = :km,
    time_unit::Symbol = :h,
    x_tick_format::AbstractString = "{:.0f}",
    y_tick_format::AbstractString = "{:.0f}",
    z_tick_format::AbstractString = "{:.0f}",
)
```

Create contour plots of the dataset `variable` in `data`, display it and save it to `file`.

# Arguments

  - `plot_file`: File to save the plots to.

  - `data_file`: HDF5 file with PinCFlow.jl output data.

  - `fields`: Either tuples of a variable name and a temporal index (for 2D data), or tuples of a variable name, three fractions which define the relative positions of the ``\\hat{y}``-``\\hat{z}``, ``\\hat{x}``-``\\hat{z}``, and ``\\hat{x}``-``\\hat{y}`` planes, and a temporal index (for 3D data).

# Keywords

  - `colormap_name`: Colormap of choice.

  - `color_tick_format`: Format string for the ticks of the colorbars.

  - `number`: Number of contour levels.

  - `space_unit`: Unit used for the coordinates. Must be `:km` or `:m`.

  - `time_unit`: Unit used for the time. Must be `:d`, `:h`, `:min` or `:s`.

  - `x_tick_format`: Format string for the ticks of the ``x``-axis.

  - `y_tick_format`: Format string for the ticks of the ``y``-axis.

  - `z_tick_format`: Format string for the ticks of the ``z``-axis.
"""
function plot_output end
